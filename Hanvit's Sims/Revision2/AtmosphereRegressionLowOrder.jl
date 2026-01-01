using LinearAlgebra
using NonlinearSolve
using ForwardDiff
using Statistics


using SatelliteToolboxAtmosphericModels
using JLD2
using GLMakie


# begin/block just to isolate local variables
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


function poly7(x, α)
    return α[1] * x + α[2] * x^2 + α[3] * x^3 + α[4] * x^4 + α[5] * x^5 + α[6] * x^6 + α[7] * x^7
end


function poly7_obj(α, data)
    t, y = data
    ŷ = poly7.(t, (α, ))
    r = ŷ - y
    return r
end


function regress_poly7(samples, p′_log)
    α_i::Vector{BigFloat} = BigFloat.(zeros(7))
    my_problem = NonlinearLeastSquaresProblem(NonlinearFunction(poly7_obj), α_i, (BigFloat.(samples), BigFloat.(p′_log)))
    my_solution = solve(my_problem)
    return my_solution.u
end


function low_poly(x, β)
    return β[1] * x + β[2] * x^2 + β[3] * x^3 + β[4] * x^4 + β[5] * x^5
end


function low_obj(β, data)
    t, y = data
    ŷ = low_poly.(t, (β, ))
    r = ŷ - y
    return r
end


function low_regression(samples, p′_log)
    β_i::Vector{BigFloat} = BigFloat.(zeros(5))
    my_problem = NonlinearLeastSquaresProblem(NonlinearFunction(low_obj), β_i, (BigFloat.(samples), BigFloat.(p′_log)))
    my_solution = solve(my_problem)
    return my_solution.u
end


function poly7_pressure(h::BigFloat, coeffs::Vector{BigFloat})
    local p::BigFloat
    if h > 5.0 * 10.0^4
        p = coeffs[12] * exp(coeffs[11] * (h - 5.0 * 10.0^4))
    elseif h < -50.0
        p = coeffs[10] * exp(coeffs[9] * (h - (-50.0)))
    else
        p = coeffs[8] * exp(poly7(h, coeffs[1:7]))
    end
    return p
end
function poly7_pressure(h::Float64, coeffs::Vector{Float64})
    local p::Float64
    if h > 5.0 * 10.0^4
        p = coeffs[12] * exp(coeffs[11] * (h - 5.0 * 10.0^4))
    elseif h < -50.0
        p = coeffs[10] * exp(coeffs[9] * (h - (-50.0)))
    else
        p = coeffs[8] * exp(poly7(h, coeffs[1:7]))
    end
    return p
end
function poly7_pressure(h::Float32, coeffs::Vector{Float32})
    local p::Float32
    if h > Float32(5.0 * 10.0^4)
        p = coeffs[12] * exp(coeffs[11] * (h - Float32(5.0 * 10.0^4)))
    elseif h < Float32(-50.0)
        p = coeffs[10] * exp(coeffs[9] * (h - Float32(-50.0)))
    else
        p = coeffs[8] * exp(poly7(h, coeffs[1:7]))
    end
    return p
end


function get_R_squared(p, h, coeffs)
    return 1.0 - (
        sum(abs2.(p - poly7_pressure.(Float32.(h), (coeffs, ))), dims = 1)[1] /
        sum(abs2.(p .- mean(p)), dims = 1)[1]
    )
end


function lowalt_pressure(h::BigFloat, coeffs::Vector{BigFloat}, base, top)
    local p::BigFloat
    if h > top
        p = coeffs[end] * exp(coeffs[end - 1] * (h - top))
    elseif h < base
        p = coeffs[end - 2] * exp(coeffs[end - 3] * (h - base))
    else
        p = coeffs[end - 4] * exp(low_poly(h, coeffs[1:end - 5]))
    end
    return p
end
function lowalt_pressure(h::Float32, coeffs::Vector{Float32}, base, top)
    local p::Float32
    if h > Float32(top)
        p = coeffs[10] * exp(coeffs[9] * (h - Float32(top)))
    elseif h < Float32(base)
        p = coeffs[8] * exp(coeffs[7] * (h - Float32(base)))
    else
        p = coeffs[6] * exp(low_poly(h, coeffs[1:5]))
    end
    return p
end


function get_R_squared_low(p, h, coeffs, base, top)
    return 1.0 - (
        sum(abs2.(p - lowalt_pressure.(Float32.(h), (coeffs, ), (base, ), (top, ))), dims = 1)[1] /
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
polyfit::Vector{BigFloat} = regress_poly7(h_samples, p′_log)
linfit_base::BigFloat = ForwardDiff.derivative(log_atm_1d, h_base)
linfit_top::BigFloat  = ForwardDiff.derivative(log_atm_1d, h_top)
p_base::BigFloat = p0 * exp(poly7(h_base, polyfit))
p_top::BigFloat  = p0 * exp(poly7(h_top, polyfit))
poly7_coeffs::Vector{BigFloat} = [
    polyfit;
    p0;
    linfit_base;
    p_base;
    linfit_top;
    p_top
]
poly7_f32_coeffs = Float32.(poly7_coeffs)
R_sq = get_R_squared(atm_data[:, 4], h_samples, poly7_f32_coeffs)


h_base_low::BigFloat = -1000.0
h_top_low::BigFloat  = 12.0 * 10.0^3
h_samples_low::Vector{BigFloat} = collect(range(h_base_low, h_top_low, length = 2^9))
atm_data_low::Matrix{BigFloat}  = collect(reduce(hcat, get_atmosphere_data.((35.34788, ), (-117.80683, ), h_samples_low))')
p_log_low::Vector{BigFloat} = log.(atm_data_low[:, 4])
p0_low::BigFloat = get_pressure(35.34788, -117.80683, BigFloat(0.0))
p0_log_low::BigFloat = log(p0_low)
p′_log_low::Vector{BigFloat} = p_log_low .- p0_log_low
polyfit_low::Vector{BigFloat} = low_regression(h_samples_low, p′_log_low)
linfit_base_low::BigFloat = ForwardDiff.derivative(log_atm_1d, h_base_low)
linfit_top_low::BigFloat  = ForwardDiff.derivative(log_atm_1d, h_top_low)
p_base_low::BigFloat = p0_low * exp(low_poly(h_base_low, polyfit_low))
p_top_low::BigFloat  = p0_low * exp(low_poly(h_top_low, polyfit_low))
lowalt_coeffs::Vector{BigFloat} = [
    polyfit_low;
    p0_low;
    linfit_base_low;
    p_base_low;
    linfit_top_low;
    p_top_low
]
lowalt_f32_coeffs = Float32.(lowalt_coeffs)
R_sq_2 = get_R_squared_low(atm_data_low[:, 4], h_samples_low, lowalt_f32_coeffs, h_base_low, h_top_low)


save(normpath(joinpath((@__FILE__), raw"..\atm_dat.jld2")), Dict(
    "oa_bf_fit" => poly7_coeffs,
    "oa_f32_fit" => poly7_f32_coeffs,
    "lowalt_bf_fit" => lowalt_coeffs,
    "lowalt_f32_fit" => lowalt_f32_coeffs,
    "fit_R2" => [R_sq; R_sq_2]
))


#=
fig = Figure()
ax = Axis(fig[1, 1], title = "Test Case", xlabel = "Altitude (metres)", ylabel = "Pressure (Pascal)")
plot_samples::Vector{BigFloat} = collect(range(-10.0^3, 1.5 * 10.0^4, length = 2^8))
actual_values::Vector{BigFloat} = get_pressure.(35.34788, -117.80683, plot_samples)
lines!(ax, plot_samples, actual_values, label = "Raw Data")
lines!(ax, plot_samples, lowalt_pressure.(plot_samples, (low_coeffs, ), (h_base_low, ), (h_top_low, )), label = "Fit")
Legend(fig, ax)
fig


fig = Figure()
ax = Axis(fig[1, 1], title = "Test Case", xlabel = "altitude (metres)", ylabel = "log pressure")
plot_samples::Vector{BigFloat} = collect(range(-10.0^3, 1.5 * 10.0^4, length = 2^8))
actual_values::Vector{BigFloat} = log_atm_1d.(plot_samples)
lines!(ax, plot_samples, actual_values, label = "Raw Data")
#lines!(ax, plot_samples, poly7_pressure.(plot_samples, (coeffs, )), label = "Fit")
Legend(fig, ax)
fig


fig = Figure()
ax = Axis(fig[1, 1], title = "Test Case", xlabel = "Altitude (metres)", ylabel = "Pressure (Pascal)")
plot_samples::Vector{BigFloat} = collect(range(-10.0^3, 10.0 * 10.0^4, length = 2^8))
actual_values::Vector{BigFloat} = get_pressure.(35.34788, -117.80683, plot_samples)
lines!(ax, plot_samples, actual_values, label = "Raw Data")
lines!(ax, plot_samples, poly7_pressure.(plot_samples, (coeffs, )), label = "Fit")
Legend(fig, ax)
fig
=#






