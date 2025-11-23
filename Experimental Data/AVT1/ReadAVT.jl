using LinearAlgebra
using NonlinearSolve

using CSV
using DataFrames
using GLMakie


include("InverseBarometerMeasurement.jl")


function fit_constant(α, data)
    x, y = data
    r = y .- α[1]
    return r
end


function fit_linear(α, data)
    x, y = data
    r = (α[1] * x .+ α[2]) - y
    return r
end


function get_constant_fit(x, y)
    my_problem = NonlinearLeastSquaresProblem(NonlinearFunction(fit_constant), [0.0], (x, y))
    my_solution = solve(my_problem)
    μ = (my_solution.u)[1]
    return μ
end


function get_linear_fit(x, y)
    my_problem = NonlinearLeastSquaresProblem(NonlinearFunction(fit_linear), [0.0; 0.0], (x, y))
    my_solution = solve(my_problem)
    lin_fit = my_solution.u
    return lin_fit
end


function get_mean(p, x)
    return dot(p, x)
end


function get_covariance(p, x, μ_x, y, μ_y)
    return get_mean(p, (x .- μ_x) .* (y .- μ_y))
end


function get_covariance_vec3(p, v, μ_v)
    xx = get_covariance(p, v[1, :], μ_v[1], v[1, :], μ_v[1])
    yy = get_covariance(p, v[2, :], μ_v[2], v[2, :], μ_v[2])
    zz = get_covariance(p, v[3, :], μ_v[3], v[3, :], μ_v[3])
    xy = get_covariance(p, v[1, :], μ_v[1], v[2, :], μ_v[2])
    xz = get_covariance(p, v[1, :], μ_v[1], v[3, :], μ_v[3])
    yz = get_covariance(p, v[2, :], μ_v[2], v[3, :], μ_v[3])
    return [
        xx xy xz;
        xy yy yz;
        xz yz zz
    ]
end


function get_linear_cd(x, y, μ, lin_fit)
    R2 = 1.0 - sum(abs2.(y - (lin_fit[1] * x .+ lin_fit[2])), dims = 1)[1] * sum(abs2.(y .- μ), dims = 1)[1]^-1
    return R2
end


csv = CSV.read(normpath(joinpath((@__FILE__), raw"..\AVTTest2.csv")), DataFrame)
t = collect(csv[:, 1])
Δt = t[2:end] - t[1:end - 1]
τ = Δt * sum(Δt, dims = 1)[1]^-1


a = [transpose(collect(csv[:, 2])); transpose(collect(csv[:, 3])); transpose(collect(csv[:, 4]))]
μ_a = [get_mean(τ, a[1, 2:end]); get_mean(τ, a[2, 2:end]); get_mean(τ, a[3, 2:end])]
cov_a = get_covariance_vec3(τ, a[:, 2:end], μ_a)
#=
Δa = (a[:, 2:end] - a[:, 1:end - 1])
μ_Δa = [get_mean(τ, Δa[1, :]); get_mean(τ, Δa[2, :]); get_mean(τ, Δa[3, :])]
cov_Δa = get_covariance_vec3(τ, Δa, μ_Δa)
=#
# rates come in millidegrees per second, and in the order pitch yaw roll
# millidegrees to radians is 10^-3 * π/180
# we want yaw pitch roll
ω = π * 180000.0^-1 * [transpose(collect(csv[:, 6])); transpose(collect(csv[:, 5])); transpose(collect(csv[:, 7]))]
μ_ω = [get_mean(τ, ω[1, 2:end]); get_mean(τ, ω[2, 2:end]); get_mean(τ, ω[3, 2:end])]
cov_ω = get_covariance_vec3(τ, ω[:, 2:end], μ_ω)
#=
Δω = (ω[:, 2:end] - ω[:, 1:end - 1])
μ_Δω = [get_mean(τ, Δω[1, :]); get_mean(τ, Δω[2, :]); get_mean(τ, Δω[3, :])]
cov_Δω = get_covariance_vec3(τ, Δω, μ_Δω)
=#


m = [transpose(collect(csv[:, 8])); transpose(collect(csv[:, 9])); transpose(collect(csv[:, 10]))]
m_norm = sqrt.(sum(abs2.(m), dims = 1))[1, :]
m̂ = m .* transpose(m_norm.^-1)
μ_m = [get_mean(τ, m[1, 2:end]); get_mean(τ, m[2, 2:end]); get_mean(τ, m[3, 2:end])]
cov_m = get_covariance_vec3(τ, m[:, 2:end], μ_m)


T = collect(csv[:, 11])
p = collect(csv[:, 12]) * 1000.0
μ_p = get_mean(τ, p[2:end])
var_p = get_covariance(τ, p[2:end], μ_p, p[2:end], μ_p)
h_baro = p2alt.(Float32.(p))
μ_h_baro = get_mean(τ, h_baro[2:end])
var_h_baro = get_covariance(τ, h_baro[2:end], μ_h_baro, h_baro[2:end], μ_h_baro)


GPS_csv = CSV.read(normpath(joinpath((@__FILE__), raw"..\GPSdata.csv")), DataFrame)
t_GPS = collect(GPS_csv[:, 7])
Δt_GPS = t_GPS[2:end] - t_GPS[1:end - 1]
τ_GPS = Δt_GPS * sum(Δt_GPS, dims = 1)[1]^-1


φ = collect(GPS_csv[:, 1])
μ_φ = get_mean(τ_GPS, φ[2:end])
var_φ = get_covariance(τ_GPS, φ[2:end], μ_φ, φ[2:end], μ_φ)
λ = collect(GPS_csv[:, 2])
μ_λ = get_mean(τ_GPS, λ[2:end])
var_λ = get_covariance(τ_GPS, λ[2:end], μ_λ, λ[2:end], μ_λ)
h = collect(GPS_csv[:, 3])
μ_h = get_mean(τ_GPS, h[2:end])
var_h = get_covariance(τ_GPS, h[2:end], μ_h, h[2:end], μ_h)


v = [transpose(collect(GPS_csv[:, 4])); transpose(collect(GPS_csv[:, 5])); transpose(collect(GPS_csv[:, 6]))]
μ_v = [get_mean(τ_GPS, v[1, 2:end]); get_mean(τ_GPS, v[2, 2:end]); get_mean(τ_GPS, v[3, 2:end])]
cov_v = get_covariance_vec3(τ_GPS, v[:, 2:end], μ_v)


#=
println("accelerometer covariance")
println(cov_a)
println("gyro covariance")
println(cov_ω)
println("magnetometer R, no hard iron")
println(cov_m)
println("barometer variance")
println(var_p)
println("GPS latitude variance")
println(var_φ)
println("GPS longitude variance")
println(var_λ)
println("GPS altitude variance")
println(var_h)
println("GPS velocity covariance")
println(cov_v)


fig = Figure()
ax = Axis(fig[1, 1], title = "Accelerometer", xlabel = "Time (seconds)", ylabel = "Acceleration (m * s^-2)")
lines!(ax, t, a[1, :], label = "x")
lines!(ax, t, a[2, :], label = "y")
lines!(ax, t, a[3, :], label = "z")
Legend(fig, ax)
fig


fig = Figure()
ax = Axis(fig[1, 1], title = "Gyroscope", xlabel = "Time (seconds)", ylabel = "Angular Velocity (rad * s^-1)")
lines!(ax, t, ω[1, :], label = "Yaw")
lines!(ax, t, ω[2, :], label = "Pitch")
lines!(ax, t, ω[3, :], label = "Roll")
Legend(fig, ax)
fig


fig = Figure()
ax = Axis(fig[1, 1], title = "Magnetometer", xlabel = "Time (seconds)", ylabel = "Magnetic Field Component (Gauss)")
lines!(ax, t, m[1, :], label = "x")
lines!(ax, t, m[2, :], label = "y")
lines!(ax, t, m[3, :], label = "z")
Legend(fig, ax)
fig


fig = Figure()
ax = Axis(fig[1, 1], title = "Barometer", xlabel = "Time (seconds)", ylabel = "Pressure (Pascal)")
ylims!(ax, 0.0, 100000.0)
lines!(ax, t, p, label = "Pressure")
Legend(fig, ax)
fig


fig = Figure()
ax = Axis(fig[1, 1], title = "GPS Velocity", xlabel = "Time (seconds)", ylabel = "Velocity (m * s^-1)")
lines!(ax, t_GPS, v[1, :], label = "v_x")
lines!(ax, t_GPS, v[2, :], label = "v_y")
lines!(ax, t_GPS, v[3, :], label = "v_z")
Legend(fig, ax)
fig


fig = Figure()
ax = Axis(fig[1, 1], title = "GPS Altitude", xlabel = "Time (seconds)", ylabel = "Altitude (m)")
lines!(ax, t_GPS, h, label = "GPS-measured geodetic altitude")
Legend(fig, ax)
fig


fig = Figure()
ax = Axis(fig[1, 1], title = "Altitude Measurement", xlabel = "Normalized Time", ylabel = "Altitude (m)")
lines!(ax, t, h_baro .- 75, label = "Barometer-measured geodetic altitude")
#lines!(ax, t / t[end], T, label = "Temperature")
lines!(ax, t_GPS, h, label = "GPS-measured geodetic altitude")
ylims!(ax, 0.0, 350.0)
Legend(fig, ax)
fig


fig = Figure()
ax = Axis(fig[1, 1], title = "Predicted Altitude", xlabel = "Altitude (metres)", ylabel = "Pressure (Pascals)")
scatter!(ax, h_baro, p, label = "calculated curve")
scatter!(ax, h_baro, get_pressure.(h_baro, (Float32.([5.0122185; -4.9929004e-5; -5.415637e-10; -3.837231e-14; 2.55155e-18; -5.321706e-23; 4.813401e-28; -1.6294356e-33]), )), label = "reference fit")
Legend(fig, ax)
fig
=#


