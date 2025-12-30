using LinearAlgebra
using Statistics


using JLD2
using GLMakie


function compensate_non_inertial(a::Vector{Float64}, φ::Float64, h::Float64)
    g = 9.780327 * (1.0 + 5.3024 * 10.0^-3 * sind(φ)^2 - 5.8 * 10.0^-6 * sind(2.0φ)^2) - (3.0877 * 10.0^-6 - 4.4 * 10.0^-9 * sind(φ)^2) * h + 7.2 * 10.0^-14 * h^2
    return a + [0.0; 0.0; g]
end
function compensate_non_inertial(φ::Float64, h::Float64)
    g = 9.780327 * (1.0 + 5.3024 * 10.0^-3 * sind(φ)^2 - 5.8 * 10.0^-6 * sind(2.0φ)^2) - (3.0877 * 10.0^-6 - 4.4 * 10.0^-9 * sind(φ)^2) * h + 7.2 * 10.0^-14 * h^2
    return [0.0; 0.0; g]
end


function get_eigenvalues_only(A::Matrix{Float64})
    return eigen(A).values
end


function get_mean_as_vec(A::Matrix{Float64}, n)
    return mean(A, dims = n)[:]
end


data = load(normpath(joinpath((@__FILE__), raw"..\FT1dat.jld2")))


t = data["t"]
τ_GPS = data["t_GPS"]
φ = data["lat"]
λ = data["lon"]
h = data["alt"]
a_1 = data["a_1"]
a_2 = data["a_2"]
a_3 = data["a_3"]
ω_1 = data["ω_1"] * 10.0^-3 * π / 180.0
ω_2 = data["ω_2"] * 10.0^-3 * π / 180.0
ω_3 = data["ω_3"] * 10.0^-3 * π / 180.0
ps_a = data["ps_a"]
last_time_before_30 = findlast(t .<= 30.0)
last_GPS_before_30 = findlast(τ_GPS .<= 30.0)


ā_1_before_30 = mean(a_1[:, 1:last_time_before_30], dims = 2)[:]
ā_2_before_30 = mean(a_2[:, 1:last_time_before_30], dims = 2)[:]
ā_3_before_30 = mean(a_3[:, 1:last_time_before_30], dims = 2)[:]
ω̄_1_before_30 = mean(ω_1[:, 1:last_time_before_30], dims = 2)[:]
ω̄_2_before_30 = mean(ω_2[:, 1:last_time_before_30], dims = 2)[:]
ω̄_3_before_30 = mean(ω_3[:, 1:last_time_before_30], dims = 2)[:]
φ̄_before_30 = mean(φ[1:last_GPS_before_30])
λ̄_before_30 = mean(λ[1:last_GPS_before_30])
h̄_before_30 = mean(h[1:last_GPS_before_30])
a_g = compensate_non_inertial(φ̄_before_30, h̄_before_30)
ab_1 = ā_1_before_30 + a_g
ab_2 = ā_2_before_30 + a_g
ab_3 = ā_3_before_30 + a_g
a′_1 = a_1 .- ab_1
a′_2 = a_2 .- ab_2
a′_3 = a_3 .- ab_3
ω′_1 = ω_1 .- ω̄_1_before_30
ω′_2 = ω_2 .- ω̄_2_before_30
ω′_3 = ω_3 .- ω̄_3_before_30
win = 100
as_1 = [a′_1[:, i:(i + win - 1)] for i = 1:(length(t) - win + 1)]
as_2 = [a′_2[:, i:(i + win - 1)] for i = 1:(length(t) - win + 1)]
as_3 = [a′_3[:, i:(i + win - 1)] for i = 1:(length(t) - win + 1)]
ās_1 = get_mean_as_vec.(as_1, 2)
ās_2 = get_mean_as_vec.(as_2, 2)
ās_3 = get_mean_as_vec.(as_3, 2)
as_1_cov = cov.(as_1, dims = 2)
as_2_cov = cov.(as_2, dims = 2)
as_3_cov = cov.(as_3, dims = 2)
as_1_ev = get_eigenvalues_only.(as_1_cov)
as_2_ev = get_eigenvalues_only.(as_2_cov)
as_3_ev = get_eigenvalues_only.(as_3_cov)
ωs_1 = [ω′_1[:, i:(i + win - 1)] for i = 1:(length(t) - win + 1)]
ωs_2 = [ω′_2[:, i:(i + win - 1)] for i = 1:(length(t) - win + 1)]
ωs_3 = [ω′_3[:, i:(i + win - 1)] for i = 1:(length(t) - win + 1)]
ω̄s_1 = get_mean_as_vec.(ωs_1, 2)
ω̄s_2 = get_mean_as_vec.(ωs_2, 2)
ω̄s_3 = get_mean_as_vec.(ωs_3, 2)
ωs_1_cov = cov.(ωs_1, dims = 2)
ωs_2_cov = cov.(ωs_2, dims = 2)
ωs_3_cov = cov.(ωs_3, dims = 2)
ωs_1_ev = get_eigenvalues_only.(ωs_1_cov)
ωs_2_ev = get_eigenvalues_only.(ωs_2_cov)
ωs_3_ev = get_eigenvalues_only.(ωs_3_cov)


fig = Figure()
ax = Axis(fig[1, 1])

#=
fig = Figure()
ax = Axis(fig[1, 1], title = "Flight Test 1 Down Angular Velocity", xlabel = "time (seconds)", ylabel = "angular rate (rad * s^-1)")
lines!(ax, t, ω′_1[3, :], label = "Microcontroller 1")
lines!(ax, t, ω′_2[3, :], label = "Microcontroller 2")
lines!(ax, t, ω′_3[3, :], label = "Microcontroller 3")
axislegend(ax, position = :lb)
fig
save(normpath(joinpath((@__FILE__), raw"..\bias_correction.png")), fig)
=#

#=
fig = Figure()
ax = Axis(fig[1, 1], title = "Flight Test 1 Angular Rate", xlabel = "mean angular rate (rad * s^-1)", ylabel = "gyro noise variance (rad^2 * s^-2)")
ylims!(ax, 0.0, 0.2)
scatter!(ax, norm.(ω̄s_1), sum.(as_1_ev), markersize = 1, label = "Microcontroller 1")
scatter!(ax, norm.(ω̄s_2), sum.(as_2_ev), markersize = 1, label = "Microcontroller 2")
scatter!(ax, norm.(ω̄s_3), sum.(as_3_ev), markersize = 1, label = "Microcontroller 3")
axislegend(ax, position = :rt)
fig
save(normpath(joinpath((@__FILE__), raw"..\accelcross.png")), fig)
=#

#=
fig = Figure()
ax = Axis(fig[1, 1], title = "Flight Test 1 Acceleration", xlabel = "mean acceleration (m * s^-2)", ylabel = "accelerometer noise variance (m^2 * s^-4)")
ylims!(ax, 0.0, 0.2)
scatter!(ax, norm.(ās_1), sum.(as_1_ev), markersize = 1, label = "Microcontroller 1")
scatter!(ax, norm.(ās_2), sum.(as_2_ev), markersize = 1, label = "Microcontroller 2")
scatter!(ax, norm.(ās_3), sum.(as_3_ev), markersize = 1, label = "Microcontroller 3")
axislegend(ax, position = :rt)
fig
save(normpath(joinpath((@__FILE__), raw"..\Acceleration.png")), fig)
=#


#=
p_1 = data["p_1"]
p_2 = data["p_2"]
p_3 = data["p_3"]
p_a = data["p_a"]


fig = Figure()
ax = Axis(fig[1, 1], title = "Flight Test 1 Pressures", xlabel = "time (seconds)", ylabel = "pressure (Pascals)")
ylims!(ax, 0.0, 1.3 * 10.0^5)
#lines!(ax, t, p_1, label = "Microcontroller 1")
lines!(ax, t, p_mc2, label = "Microcontroller 2")
#lines!(ax, t, p_3, label = "Microcontroller 3")
lines!(ax, t, p_ahrs, label = "AHRS")
axislegend(ax, position = :lb)
fig
save(normpath(joinpath((@__FILE__), raw"..\Pressures.png")), fig)
=#








