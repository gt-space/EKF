using LinearAlgebra
using Statistics


using JLD2
using GLMakie
using CSV
using DataFrames


data = load(normpath(joinpath((@__FILE__), raw"..\FT1dat.jld2")))
τ_GPS = data["t_GPS"]
φ = data["lat"]
λ = data["lon"]
h = data["alt"]
v_true = data["vel"]


dat1 = CSV.read(normpath(joinpath((@__FILE__), raw"..\mc0WithoutGPS.csv")), DataFrame)
dat1_gps = CSV.read(normpath(joinpath((@__FILE__), raw"..\mc0withGPS.csv")), DataFrame)
dat1_dr = CSV.read(normpath(joinpath((@__FILE__), raw"..\mc0DeadReckoning.csv")), DataFrame)
dat2 = CSV.read(normpath(joinpath((@__FILE__), raw"..\mc1WithoutGPS.csv")), DataFrame)
dat2_gps = CSV.read(normpath(joinpath((@__FILE__), raw"..\mc1withGPS.csv")), DataFrame)
dat2_dr = CSV.read(normpath(joinpath((@__FILE__), raw"..\mc1DeadReckoning.csv")), DataFrame)
dat3 = CSV.read(normpath(joinpath((@__FILE__), raw"..\mc2WithoutGPS.csv")), DataFrame)
dat3_gps = CSV.read(normpath(joinpath((@__FILE__), raw"..\mc2withGPS.csv")), DataFrame)
dat3_dr = CSV.read(normpath(joinpath((@__FILE__), raw"..\mc2DeadReckoning.csv")), DataFrame)
t = Float64.(collect(dat1[:, 1])) .- 1765126645.34027
alt = [
    Float64.(collect(dat1[:, 8]))';
    Float64.(collect(dat2[:, 8]))';
    Float64.(collect(dat3[:, 8]))'
]
alt_gps = [
    Float64.(collect(dat1_gps[:, 8]))';
    Float64.(collect(dat2_gps[:, 8]))';
    Float64.(collect(dat3_gps[:, 8]))'
]
alt_dr = [
    Float64.(collect(dat1_dr[:, 8]))';
    Float64.(collect(dat2_dr[:, 8]))';
    Float64.(collect(dat3_dr[:, 8]))'
]
v_EKF = [
    Float64.(collect(dat1[:, 11]))';
    Float64.(collect(dat2[:, 11]))';
    Float64.(collect(dat3[:, 11]))'
]
v_gps = [
    Float64.(collect(dat1_gps[:, 11]))';
    Float64.(collect(dat2_gps[:, 11]))';
    Float64.(collect(dat3_gps[:, 11]))'
]
v_dr = [
    Float64.(collect(dat1_dr[:, 11]))';
    Float64.(collect(dat2_dr[:, 11]))';
    Float64.(collect(dat3_dr[:, 11]))'
]

φ_EKF = [
    Float64.(collect(dat1[:, 6]))';
    Float64.(collect(dat2[:, 6]))';
    Float64.(collect(dat3[:, 6]))';
]
λ_EKF = [
    Float64.(collect(dat1[:, 7]))';
    Float64.(collect(dat2[:, 7]))';
    Float64.(collect(dat3[:, 7]))';
]


t_cutoff = 310.0

fig = Figure()
ax = Axis(fig[1, 1:2], title = "EKF Latitude", xlabel = "time (seconds)", ylabel = "latitude (degrees)")
lines!(ax, t[1:findlast(t .<= t_cutoff)], φ_EKF[1, 1:findlast(t .<= t_cutoff)], label = "Microcontroller 1")
lines!(ax, t[1:findlast(t .<= t_cutoff)], φ_EKF[2, 1:findlast(t .<= t_cutoff)], label = "Microcontroller 2")
lines!(ax, t[1:findlast(t .<= t_cutoff)], φ_EKF[3, 1:findlast(t .<= t_cutoff)], label = "Microcontroller 3")
lines!(ax, τ_GPS[1:findlast(τ_GPS .<= t_cutoff)], φ[1:findlast(τ_GPS .<= t_cutoff)], label = "Reference (FC GPS)")
axislegend(ax, position = :lb)

ax2 = Axis(fig[2, 1:2], title = "EKF Longitude", xlabel = "time (seconds)", ylabel = "latitude (degrees)")
lines!(ax2, t[1:findlast(t .<= t_cutoff)], λ_EKF[1, 1:findlast(t .<= t_cutoff)], label = "Microcontroller 1")
lines!(ax2, t[1:findlast(t .<= t_cutoff)], λ_EKF[2, 1:findlast(t .<= t_cutoff)], label = "Microcontroller 2")
lines!(ax2, t[1:findlast(t .<= t_cutoff)], λ_EKF[3, 1:findlast(t .<= t_cutoff)], label = "Microcontroller 3")
lines!(ax2, τ_GPS[1:findlast(τ_GPS .<= t_cutoff)], λ[1:findlast(τ_GPS .<= t_cutoff)], label = "Reference (FC GPS)")
axislegend(ax2, position = :lb)
fig

#=
fig = Figure()
ax = Axis(fig[1, 1], title = "EKF With GPS", xlabel = "time (seconds)", ylabel = "down velocity (metres * seconds^-1)")
lines!(ax, t, v_EKF[1, :], label = "Microcontroller 1")
lines!(ax, t, v_EKF[2, :], label = "Microcontroller 2")
lines!(ax, t, v_EKF[3, :], label = "Microcontroller 3")
lines!(ax, τ_GPS, v_true[3, :], label = "Reference (FC GPS Down Velocity)")
axislegend(ax, position = :lt)
fig
save(normpath(joinpath((@__FILE__), raw"..\ekf_vel_no_gps_no_fault.png")), fig)
fig
=#

