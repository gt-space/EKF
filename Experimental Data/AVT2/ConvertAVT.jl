using XLSX
using JLD2


xf = XLSX.readxlsx(normpath(joinpath((@__FILE__), raw"..\ShortAVT2.xlsx")))
sh = xf["AVT"]


timestamps = Float64.(collect(sh[2:end, 1]))[:]


φ = Float64.(collect(sh[2:end, 5]))[:]
λ = Float64.(collect(sh[2:end, 6]))[:]
h = Float64.(collect(sh[2:end, 2]))[:]
GPS_valid = (collect(sh[2:end, 8]))[:] .== "true"
v_x = Float64.(collect(sh[2:end, 7]))[:]
v_y = Float64.(collect(sh[2:end, 4]))[:]
v_z = Float64.(collect(sh[2:end, 3]))[:]


ω_x = Float64.(collect(sh[2:end, 9]))[:]
ω_y = Float64.(collect(sh[2:end, 10]))[:]
ω_z = Float64.(collect(sh[2:end, 11]))[:]
a_x = Float64.(collect(sh[2:end, 12]))[:]
a_y = Float64.(collect(sh[2:end, 13]))[:]
a_z = Float64.(collect(sh[2:end, 14]))[:]


m_x = Float64.(collect(sh[2:end, 15]))[:]
m_y = Float64.(collect(sh[2:end, 16]))[:]
m_z = Float64.(collect(sh[2:end, 17]))[:]


p = Float64.(collect(sh[2:end, 18]))[:]
T = Float64.(collect(sh[2:end, 19]))[:]


save(normpath(joinpath((@__FILE__), raw"..\AVT_data2.jld2")), Dict(
    "times" => timestamps,
    "lat" => φ,
    "lon" => λ,
    "alt" => h,
    "valid" => GPS_valid,
    "v_x" => v_x,
    "v_y" => v_y,
    "v_z" => v_z,
    "ω_x" => ω_x,
    "ω_y" => ω_y,
    "ω_z" => ω_z,
    "a_x" => a_x,
    "a_y" => a_y,
    "a_z" => a_z,
    "m_x" => m_x,
    "m_y" => m_y,
    "m_z" => m_z,
    "pressure" => p,
    "temp" => T   
))




