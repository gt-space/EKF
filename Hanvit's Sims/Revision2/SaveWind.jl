using JLD2
using CSV
using DataFrames


wind_csv = CSV.read(normpath(joinpath((@__FILE__), raw"..\test_wind_data.csv")), DataFrame)
wind_sp = reverse(collect(wind_csv[:, 1]))
wind_v  = transpose([reverse(collect(wind_csv[:, 2])) reverse(collect(wind_csv[:, 3])) reverse(collect(wind_csv[:, 4]))])[:, :]


save(normpath(joinpath((@__FILE__), raw"..\test_wind_data.jld2")), Dict(
    "sample_points" => wind_sp,
    "velocity_vectors" => wind_v
))





