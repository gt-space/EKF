using JLD2
using CSV
using DataFrames


CAon_csv = CSV.read(normpath(joinpath((@__FILE__), raw"..\CAon.csv")), DataFrame)
CAoff_csv = CSV.read(normpath(joinpath((@__FILE__), raw"..\CAoff.csv")), DataFrame)
CN_csv = CSV.read(normpath(joinpath((@__FILE__), raw"..\CN.csv")), DataFrame)
CM_csv = CSV.read(normpath(joinpath((@__FILE__), raw"..\CM.csv")), DataFrame)


CAon::Matrix{Float64} = zeros(450, 16)
CAoff::Matrix{Float64} = zeros(450, 16)
for i = 0:15
    CAon[:, i + 1] = collect(CAon_csv[:, i + 2])
    CAoff[:, i + 1] = collect(CAoff_csv[:, i + 2])
end
CN::Vector{Float64} = collect(CN_csv[:, 2])
CM::Vector{Float64} = collect(CM_csv[:, 2])


Atmosphere_csv = CSV.read(normpath(joinpath((@__FILE__), raw"..\Atmosphere.csv")), DataFrame)


heights::Vector{Float64} = collect(Atmosphere_csv[1, :]) * 10.0^3
densities::Vector{Float64} = collect(Atmosphere_csv[2, :])
temperatures::Vector{Float64} = collect(Atmosphere_csv[3, :])
molar_masses::Vector{Float64} = collect(Atmosphere_csv[4, :])
pressures::Vector{Float64} = densities .* (8.31446261815324 * molar_masses.^-1) .* temperatures


save(normpath(joinpath((@__FILE__), raw"..\PerformanceData.jld2")), Dict(
    "CAon" => CAon,
    "CAoff" => CAoff,
    "CN" => CN,
    "CM" => CM,
    "heights" => heights,
    "densities" => densities,
    "temperatures" => temperatures,
    "molar_masses" => molar_masses,
    "pressures" => pressures
))
# my_dict = load("PerformanceData.jld2")



