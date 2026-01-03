using LinearAlgebra
using NonlinearSolve
using ForwardDiff
using Statistics


using SatelliteToolboxAtmosphericModels
using JLD2
using GLMakie


begin
    # FAR coordinates: 35.34788, -117.80683, 625.0
    t_launch = DateTime("2025-01-10T11:00:0.0")
    φ_launch = 35.34788
    λ_launch = -117.80683
    h_launch = 625.0

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
        pressure_nrlmsise00(h)

    Return the barometric pressure of the atmosphere at the specified geodetic altitude.

    ...
    # Arguments
    - `h::Real`: Geodetic altitude (metres)

    ...
    """
    function pressure_nrlmsise00(h)
        local atm_data = AtmosphericModels.nrlmsise00(
            t_launch,
            h,
            φ_launch,
            λ_launch
        )
        local frac_species = normalized_abundances(atm_data)
        local T = atm_data.temperature
        local m = dot(molar_masses, frac_species)
        local ρ = atm_data.total_density
        return p = ρ * (8.31446261815324 * m^-1 * T)
    end


    """
        true_log_pressure(h)

    Return the natural logarithm of the barometric pressure of the atmosphere at the specified geodetic altitude.

    ...
    # Arguments
    - `h::Real`: Geodetic altitude (metres)

    ...
    """
    function true_log_pressure(h)
        return log(pressure_nrlmsise00(h))
    end


    """
        true_log_deriv_pressure(h)

    Return the logarithmic derivative with respect to geodetic altitude of the barometric pressure of the atmosphere at the specified geodetic altitude.

    ...
    # Arguments
    - `h::Real`: Geodetic altitude (metres)

    ...
    """
    function true_log_deriv_pressure(h)
        return ForwardDiff.derivative(true_log_pressure, h)
    end


    """
        true_log_deriv_pressure(h)

    Return the second derivative with respect to geodetic altitude of the logarithm of the barometric pressure of the atmosphere at the specified geodetic altitude.

    ...
    # Arguments
    - `h::Real`: Geodetic altitude (metres)

    ...
    """
    function true_log_2nd_deriv_pressure(h)
        return ForwardDiff.derivative(true_log_deriv_pressure, h)
    end
end


"""
    poly5(x, β)

Calculate the value of a fifth-order polynomial with no constant term

...
# Arguments
- `x::Real`: Polynomial argument
- `β::Vector{Real}`: Five-element vector containing the polynomial fit coefficients from 1st to 5th order

...
"""
function poly5(x, β)
    return β[1] * x + β[2] * x^2 + β[3] * x^3 + β[4] * x^4 + β[5] * x^5
end


"""
    poly5_res(β, data)

Find the residuals of the polynomial interpolant

...
# Arguments
- `β::Vector{Real}`: Five-element vector containing the polynomial fit coefficients from 1st to 5th order
- `data::Tuple{Vector{BigFloat}, Vector{BigFloat}}`: Pressure dataset and geodetic altitude samples to regress against.

...
"""
function poly5_res(β, data::Tuple{Vector{BigFloat}, Vector{BigFloat}})
    t_true, y_true = data
    ŷ = poly5.(t_true, (β, ))
    r = ŷ - y_true
    return r
end
function poly5_res(β::Vector{Float32}, data::Tuple{Vector{BigFloat}, Vector{BigFloat}})
    t_true, y_true = data
    ŷ::Vector{Float32} = poly5.(Float32.(t_true), (β, ))
    r::Vector{BigFloat} = BigFloat.(ŷ) - y_true
    return r
end


"""
    regress_poly5(data::Tuple{Vector{BigFloat}, Vector{BigFloat}})

Fit a fifth order polynomial to the provided data

...
# Arguments
- `data::Tuple{Vector{BigFloat}, Vector{BigFloat}}`: Pressure dataset and geodetic altitude samples to regress against.

...
"""
function regress_poly5(data::Tuple{Vector{BigFloat}, Vector{BigFloat}})
    β_0::Vector{BigFloat} = BigFloat.(zeros(5))
    my_problem = NonlinearLeastSquaresProblem(NonlinearFunction(poly5_res), β_0, data)
    my_solution = solve(my_problem)
    return my_solution.u
end


"""
    tsos_poly5(x::Vector{T}) where T <: Real

Find the sum of squares of the given dataset

...
# Arguments
- `x::Vector{T}`: Dataset to find the total sum-of-squares of

...
"""
function tsos_poly5(x::Vector{T}) where T <: Real
    x̄ = mean(x)
    buffer = x .- x̄
    return dot(buffer, buffer)
end


"""
    rsq_poly5(data::Tuple{Vector{T}, Vector{T}}) where T<: Real

Find the coefficient of determination of the given dataset

...
# Arguments
- `α::Vector{Real}`: Seven-element vector containing the polynomial fit coefficients from 1st to 5th order
- `data::Tuple{Vector{BigFloat}, Vector{BigFloat}}`: Pressure dataset and geodetic altitude samples to regress against.

...
"""
function rsq_poly5(data::Tuple{Vector{T}, Vector{T}}) where T<: Real
    t_true, y_true = data
    res = poly5_res(Float32.(β), data)
    return 1.0 - dot(res, res) * tsos_poly5(y_true)^-1
end


function logatm5(h)
    local logp
    if h < h_base
        logp = m_base * (h - h_base) + b_base
    elseif h > h_ceil
        logp = m_ceil * (h - h_ceil) + b_ceil
    else
        logp = poly5(h, β)
    end
    return logp
end
function logatm5(h::Float32)
    local logp::Float32
    if h < Float32(h_base)
        logp = Float32(m_base) * (h - Float32(h_base)) + Float32(b_base)
    elseif h > Float32(h_ceil)
        logp = Float32(m_ceil) * (h - h_ceil) + Float32(b_ceil)
    else
        logp = poly5(h, Float32.(β))
    end
    return logp
end


setprecision(128)
SpaceIndices.init()


h_base::BigFloat = -10.0^3
h_ceil::BigFloat = 12.0 * 10.0^3
h::Vector{BigFloat}   = collect(range(h_base, h_ceil, length = 2^9))
p::Vector{BigFloat}   = true_log_pressure.(h) .- true_log_pressure(0.0)
p′::Vector{BigFloat}  = true_log_deriv_pressure.(h)
p′′::Vector{BigFloat} = true_log_2nd_deriv_pressure.(h)
β::Vector{BigFloat}   = regress_poly5((h, p))
m_base::BigFloat      = true_log_deriv_pressure(h_base)
m_ceil::BigFloat      = true_log_deriv_pressure(h_ceil)
b_base::BigFloat      = poly5(h_base, β)
b_ceil::BigFloat      = poly5(h_ceil, β)
P_0::BigFloat         = pressure_nrlmsise00(0.0)
C_base::BigFloat      = m_base * P_0
C_ceil::BigFloat      = m_ceil * P_0


h_plot::Vector{BigFloat}   = collect(range(-1.0 * 10.0^3, 12.0 * 10.0^3, length = 2^11))
p_plot::Vector{BigFloat}   = true_log_pressure.(h_plot) .- true_log_pressure(0.0)
p_interp::Vector{BigFloat} = BigFloat.(logatm5.(Float32.(h_plot)))
rsq::BigFloat = rsq_poly5((h_plot, p_plot))
ε_log::Vector{BigFloat} = p_interp - p_plot
ε_log ./= p_plot
ε_log_rms::BigFloat = sqrt(mean(abs2.(ε_log)))
ε_log_max::BigFloat = maximum(abs.(ε_log))
p_lin::Vector{BigFloat} = pressure_nrlmsise00.(h_plot)
ε_lin::Vector{BigFloat} = P_0 * exp.(p_interp) - p_lin
ε_lin ./= p_lin
ε_lin_rms::BigFloat = sqrt(mean(abs2.(ε_lin)))
ε_lin_max::BigFloat = maximum(abs.(ε_lin))
#=
fig = Figure(size = (3000, 2000))
ax = Axis(fig[1, 1])
lines!(ax, h_plot, p_plot, label = "true log pressure")
lines!(ax, h_plot, p_interp, label = "interpolated log pressure")
axislegend(ax, position = :lb)
fig
save(normpath(joinpath((@__FILE__), raw"..\altimeterbarometerfunctiondiagnostic.png")), fig, size = (3000, 2000))
fig
=#


save(normpath(joinpath((@__FILE__), raw"..\AltimeterBarometerData.jld2")), Dict(
    "h_base" => h_base,
    "h_ceil" => h_ceil,
    "m_base" => m_base,
    "m_ceil" => m_ceil,
    "b_base" => b_base,
    "b_ceil" => b_ceil,
    "C_base" => C_base,
    "C_ceil" => C_ceil,
    "β"      => β,
    "P_0"    => P_0,
    "rsq"    => rsq,
    "rms_rel_error"   => ε_lin_rms,
    "worst_rel_error" => ε_lin_max
))


save(normpath(joinpath((@__FILE__), raw"..\AltimeterBarometerDataF32.jld2")), Dict(
    "h_base" => Float32(h_base),
    "h_ceil" => Float32(h_ceil),
    "m_base" => Float32(m_base),
    "m_ceil" => Float32(m_ceil),
    "b_base" => Float32(b_base),
    "b_ceil" => Float32(b_ceil),
    "C_base" => Float32(C_base),
    "C_ceil" => Float32(C_ceil),
    "β"      => Float32.(β),
    "P_0"    => Float32(P_0)
))

