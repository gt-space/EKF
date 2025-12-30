begin


local wind_file = load(normpath(joinpath((@__FILE__), raw"..\test_wind_data.jld2")))
local wind_sp = wind_file["sample_points"]
local wind_v = wind_file["velocity_vectors"]


function get_wind(φ::Float64, λ::Float64, h::Float64)
    local w::Vector{Float64} = lerp(h, wind_sp, wind_v)
    return w
end


end


begin


# molar masses from PubChem
local molar_masses::Vector{Float64} = [
    0.01400674;
    0.028014;
    0.0159994;
    0.031999;
    0.00100794;
    0.004002602; 
    0.039948
]


local function normalized_abundances(atm_data)
    local Σ_number_density::Float64 = (
        atm_data.N_number_density +
        atm_data.N2_number_density +
        atm_data.O_number_density + atm_data.aO_number_density +
        atm_data.O2_number_density +
        atm_data.H_number_density +
        atm_data.He_number_density + 
        atm_data.Ar_number_density
    )
    return abundances::Vector{Float64} = Σ_number_density^-1 * [
        atm_data.N_number_density;
        atm_data.N2_number_density;
        atm_data.O_number_density + atm_data.aO_number_density;
        atm_data.O2_number_density;
        atm_data.H_number_density;
        atm_data.He_number_density;
        atm_data.Ar_number_density
    ]
end


"""
    get_atmosphere_data(φ::Float64, λ::Float64, h::Float64)

Return temperature, molar mass, density, and pressure of the atmosphere at the specified LLA coordinates.

...
# Arguments
- `φ::Float64`: Geodetic latitude
- `λ::Float64`: Geodetic longitude
- `h::Float64`: Geodetic altitude

...
"""
function get_atmosphere_data(φ::Float64, λ::Float64, h::Float64)
    # FAR coordinates: 35.34788, -117.80683, 625.0
    local atm_data = AtmosphericModels.nrlmsise00(
        (DateTime("2025-01-10T11:00:0.0")),
        h,
        φ,
        λ
    )
    local P_species::Vector{Float64} = normalized_abundances(atm_data)
    local T::Float64 = atm_data.temperature
    local m::Float64 = dot(molar_masses, P_species)
    local ρ::Float64 = atm_data.total_density
    local p::Float64 = ρ * ((8.31446261815324 * m^-1) * T)
    return atm_properties::Vector{Float64} = [
        T;
        m;
        ρ;
        p
    ]
end


end

#=
h_samples = collect(range(-1000.0, 65.0 * 10.0^3, length = 10^5))
atm_mat = collect(reduce(hcat, get_atmosphere_data.((35.34788, ), (-117.80683, ), h_samples))')


using GLMakie


fig = Figure()
ax = Axis(fig[1, 1], title = "Temperature Plot", xlabel = "Altitude (metres)", ylabel = "Temperature (Kelvin)")
lines!(ax, h_samples, atm_mat[:, 1])
fig
=#

#=
h_samples = [collect(range(-10.0^2, 5.0 * 10.0^3, length = 2^11+1))[1:end - 1]; 5.0 * 10.0^3 * 10.0.^(collect(range(0.0, 1.0, length = 2^12)))]
atm_mat = collect(reduce(hcat, get_atmosphere_data.((35.34788, ), (-117.80683, ), h_samples))')
fig = Figure()
ax = Axis(fig[1, 1], title = "Pressure Plot", xlabel = "Altitude (metres)", ylabel = "Pressure (Pascal)")
lines!(ax, h_samples, log10.(atm_mat[:, 4]))
fig
=#

