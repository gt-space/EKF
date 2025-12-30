using LinearAlgebra
using NonlinearSolve
using ForwardDiff
using Statistics


using SatelliteToolboxAtmosphericModels
using JLD2
using GLMakie


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
    local Σ_number_density = (
        atm_data.N_number_density +
        atm_data.N2_number_density +
        atm_data.O_number_density + atm_data.aO_number_density +
        atm_data.O2_number_density +
        atm_data.H_number_density +
        atm_data.He_number_density + 
        atm_data.Ar_number_density
    )
    return abundances = Σ_number_density^-1 * [
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
function get_atmosphere_data(φ, λ, h)
    # FAR coordinates: 35.34788, -117.80683, 625.0
    local atm_data = AtmosphericModels.nrlmsise00(
        (DateTime("2025-01-10T11:00:0.0")),
        h,
        φ,
        λ
    )
    local P_species = normalized_abundances(atm_data)
    local T = atm_data.temperature
    local m = dot(molar_masses, P_species)
    local ρ = atm_data.total_density
    local p = ρ * ((8.31446261815324 * m^-1) * T)
    return atm_properties = [
        T;
        m;
        ρ;
        p
    ]
end


function get_pressure(φ, λ, h)
    # FAR coordinates: 35.34788, -117.80683, 625.0
    local atm_data = AtmosphericModels.nrlmsise00(
        (DateTime("2025-01-10T11:00:0.0")),
        h,
        φ,
        λ
    )
    local P_species = normalized_abundances(atm_data)
    local T = atm_data.temperature
    local m = dot(molar_masses, P_species)
    local ρ = atm_data.total_density
    return p = ρ * ((8.31446261815324 * m^-1) * T)
end


function log_atm_1d(h)
    return log(get_pressure(35.34788, -117.80683, h))
end


end


function poly8(x, α)
    return α[1] * x + α[2] * x^2 + α[3] * x^3 + α[4] * x^4 + α[5] * x^5 + α[6] * x^6 + α[7] * x^7 + α[8] * x^8
end


function poly8_obj(α, data)
    t, y = data
    ŷ = poly8.(t, (α, ))
    r = ŷ - y
    return r
end


function regress_poly8(samples, p′_log)
    β_i::Vector{BigFloat} = BigFloat.(zeros(8))
    my_problem = NonlinearLeastSquaresProblem(NonlinearFunction(poly8_obj), β_i, (BigFloat.(samples), BigFloat.(p′_log)))
    my_solution = solve(my_problem)
    return my_solution.u
end


function poly8_pressure(h::BigFloat, coeffs::Vector{BigFloat})
    local p::BigFloat
    if h > 5.0 * 10.0^4
        p = coeffs[13] * exp(coeffs[12] * (h - 5.0 * 10.0^4))
    elseif h < -50.0
        p = coeffs[11] * exp(coeffs[10] * (h - (-50.0)))
    else
        p = coeffs[9] * exp(poly8(h, coeffs[1:8]))
    end
    return p
end
function poly8_pressure(h::Float32, coeffs::Vector{Float32})
    local p::Float32
    if h > Float32(5.0 * 10.0^4)
        p = coeffs[13] * exp(coeffs[12] * (h - Float32(5.0 * 10.0^4)))
    elseif h < Float32(-50.0)
        p = coeffs[11] * exp(coeffs[10] * (h - Float32(-50.0)))
    else
        p = coeffs[9] * exp(poly8(h, coeffs[1:8]))
    end
    return p
end


function get_R_squared(p, h, coeffs)
    return 1.0 - (
        sum(abs2.(p - poly8_pressure.(Float32.(h), (coeffs, ))), dims = 1)[1] /
        sum(abs2.(p .- mean(p)), dims = 1)[1]
    )
end


setprecision(128)
SpaceIndices.init()


h_base::BigFloat = -50.0 
h_top::BigFloat  = 5.0 * 10.0^4
h_samples::Vector{BigFloat} = collect(range(h_base, h_top, length = 2^9))
atm_data::Matrix{BigFloat}  = collect(reduce(hcat, get_atmosphere_data.((35.34788, ), (-117.80683, ), h_samples))')
p_log::Vector{BigFloat} = log.(atm_data[:, 4])
p0::BigFloat = get_pressure(35.34788, -117.80683, BigFloat(0.0))
p0_log::BigFloat = log(p0)
p′_log::Vector{BigFloat} = p_log .- p0_log
polyfit::Vector{BigFloat} = regress_poly8(h_samples, p′_log)
linfit_base::BigFloat = ForwardDiff.derivative(log_atm_1d, h_base)
linfit_top::BigFloat  = ForwardDiff.derivative(log_atm_1d, h_top)
p_base::BigFloat = get_pressure(35.34788, -117.80683, h_base) 
p_top::BigFloat  = get_pressure(35.34788, -117.80683, h_top)
coeffs::Vector{BigFloat} = [
    polyfit;
    p0;
    linfit_base;
    p_base;
    linfit_top;
    p_top
]
poly8_f32_coeffs = Float32.(coeffs)
get_R_squared(atm_data[:, 4], h_samples, poly8_f32_coeffs)










