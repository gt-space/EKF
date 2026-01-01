using LinearAlgebra


using CSV
using DataFrames
using JLD2
using GLMakie


dat1 = CSV.read(normpath(joinpath((@__FILE__), raw"..\165725.csv")), DataFrame)
dat2 = CSV.read(normpath(joinpath((@__FILE__), raw"..\170854.csv")), DataFrame)
n1 = 135500
n2 = 48400
t = [Float64.(collect(dat1[:, 1][1:n1])); Float64.(collect(dat2[:, 1][1:n2]))]
println(t[1])
t .-= t[1]
GPS_valid = [collect(dat1[:, 42][1:n1]) .== true; collect(dat2[:, 42][1:n2]) .== true]
τ_GPS = t[GPS_valid]
φ_GPS = [Float64.(collect(dat1[:, 38][1:n1])); Float64.(collect(dat2[:, 38][1:n2]))][GPS_valid]
λ_GPS = [Float64.(collect(dat1[:, 39][1:n1])); Float64.(collect(dat2[:, 39][1:n2]))][GPS_valid]
h_GPS = [Float64.(collect(dat1[:, 34][1:n1])); Float64.(collect(dat2[:, 34][1:n2]))][GPS_valid]
v_GPS = [
    [Float64.(collect(dat1[:, 40][1:n1]))' Float64.(collect(dat2[:, 40][1:n2]))'][GPS_valid]';
    [Float64.(collect(dat1[:, 36][1:n1]))' Float64.(collect(dat2[:, 36][1:n2]))'][GPS_valid]';
    [Float64.(collect(dat1[:, 35][1:n1]))' Float64.(collect(dat2[:, 35][1:n2]))'][GPS_valid]'
]
p_ahrs = 1000.0 * [Float64.(collect(dat1[:, 7][1:n1])); Float64.(collect(dat2[:, 7][1:n2]))]
p_mc1 = [Float64.(collect(dat1[:, 67][1:n1])); Float64.(collect(dat2[:, 67][1:n2]))]
p_mc2 = [Float64.(collect(dat1[:, 113][1:n1])); Float64.(collect(dat2[:, 113][1:n2]))]
p_mc3 = [Float64.(collect(dat1[:, 159][1:n1])); Float64.(collect(dat2[:, 159][1:n2]))]
a_mc1 = [
    Float64.(collect(dat1[:, 58][1:n1]))' Float64.(collect(dat2[:, 58][1:n2]))';
    Float64.(collect(dat1[:, 59][1:n1]))' Float64.(collect(dat2[:, 59][1:n2]))';
    Float64.(collect(dat1[:, 60][1:n1]))' Float64.(collect(dat2[:, 60][1:n2]))'
]
a_mc2 = [
    Float64.(collect(dat1[:, 104][1:n1]))' Float64.(collect(dat2[:, 104][1:n2]))';
    Float64.(collect(dat1[:, 105][1:n1]))' Float64.(collect(dat2[:, 105][1:n2]))';
    Float64.(collect(dat1[:, 106][1:n1]))' Float64.(collect(dat2[:, 106][1:n2]))'
]
a_mc3 = [
    Float64.(collect(dat1[:, 150][1:n1]))' Float64.(collect(dat2[:, 150][1:n2]))';
    Float64.(collect(dat1[:, 151][1:n1]))' Float64.(collect(dat2[:, 151][1:n2]))';
    Float64.(collect(dat1[:, 152][1:n1]))' Float64.(collect(dat2[:, 152][1:n2]))'
]
ω_mc1 = [
    Float64.(collect(dat1[:, 49][1:n1]))' Float64.(collect(dat2[:, 49][1:n2]))';
    Float64.(collect(dat1[:, 50][1:n1]))' Float64.(collect(dat2[:, 50][1:n2]))';
    Float64.(collect(dat1[:, 51][1:n1]))' Float64.(collect(dat2[:, 51][1:n2]))'
]
ω_mc2 = [
    Float64.(collect(dat1[:, 95][1:n1]))' Float64.(collect(dat2[:, 95][1:n2]))';
    Float64.(collect(dat1[:, 96][1:n1]))' Float64.(collect(dat2[:, 96][1:n2]))';
    Float64.(collect(dat1[:, 97][1:n1]))' Float64.(collect(dat2[:, 97][1:n2]))'
]
ω_mc3 = [
    Float64.(collect(dat1[:, 141][1:n1]))' Float64.(collect(dat2[:, 141][1:n2]))';
    Float64.(collect(dat1[:, 142][1:n1]))' Float64.(collect(dat2[:, 142][1:n2]))';
    Float64.(collect(dat1[:, 143][1:n1]))' Float64.(collect(dat2[:, 143][1:n2]))'
]
win = 100
t_slices = [[t[i:i + win - 1]] for i = 1:(length(t) - win + 1)]
p_slices_ahrs = [[p_ahrs[i:i + win - 1]] for i = 1:(length(t) - win + 1)]
a_slices_mc1 = [[a_mc1[:, i:(i + win - 1)]] for i = 1:(length(t) - win + 1)]
a_slices_mc2 = [[a_mc2[:, i:(i + win - 1)]] for i = 1:(length(t) - win + 1)]
a_slices_mc3 = [[a_mc3[:, i:(i + win - 1)]] for i = 1:(length(t) - win + 1)]
ω_slices_mc1 = [[ω_mc1[:, i:(i + win - 1)]] for i = 1:(length(t) - win + 1)]
ω_slices_mc2 = [[ω_mc2[:, i:(i + win - 1)]] for i = 1:(length(t) - win + 1)]
ω_slices_mc3 = [[ω_mc3[:, i:(i + win - 1)]] for i = 1:(length(t) - win + 1)]


save(normpath(joinpath((@__FILE__), raw"..\FT1dat.jld2")), Dict(
    "t" => t,
    "GPS_valid" => GPS_valid,
    "t_GPS" => τ_GPS,
    "lat" => φ_GPS,
    "lon" => λ_GPS,
    "alt" => h_GPS,
    "vel" => v_GPS,
    "p_a" => p_ahrs,
    "p_1" => p_mc1,
    "p_2" => p_mc2,
    "p_3" => p_mc3,
    "a_1" => a_mc1,
    "a_2" => a_mc2,
    "a_3" => a_mc3,
    "ω_1" => ω_mc1,
    "ω_2" => ω_mc2,
    "ω_3" => ω_mc3,
    "ts" => t_slices,
    "ps_a" => p_slices_ahrs,
    "ωs_1" => ω_slices_mc1,
    "ωs_2" => ω_slices_mc2,
    "ωs_3" => ω_slices_mc3
))

