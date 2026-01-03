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
    poly7(x, α)

Calculate the value of a seventh-order polynomial with no constant term

...
# Arguments
- `x::Real`: Polynomial argument
- `α::Vector{Real}`: Seven-element vector containing the polynomial fit coefficients from 1st to 7th order

...
"""
function poly7(x, α)
    return α[1] * x + α[2] * x^2 + α[3] * x^3 + α[4] * x^4 + α[5] * x^5 + α[6] * x^6 + α[7] * x^7
end


"""
    poly7_res(α, data)

Find the residuals of the polynomial interpolant

...
# Arguments
- `α::Vector{Real}`: Seven-element vector containing the polynomial fit coefficients from 1st to 7th order
- `data::Tuple{Vector{BigFloat}, Vector{BigFloat}}`: Pressure dataset and geodetic altitude samples to regress against.

...
"""
function poly7_res(α, data::Tuple{Vector{BigFloat}, Vector{BigFloat}})
    t_true, y_true = data
    ŷ = poly7.(t_true, (α, ))
    r = ŷ - y_true
    return r
end
function poly7_res(α::Vector{Float32}, data::Tuple{Vector{BigFloat}, Vector{BigFloat}})
    t_true, y_true = data
    ŷ::Vector{Float32} = poly7.(Float32.(t_true), (α, ))
    r::Vector{BigFloat} = BigFloat.(ŷ) - y_true
    return r
end


"""
    regress_poly7(data::Tuple{Vector{BigFloat}, Vector{BigFloat}})

Fit a seventh order polynomial to the provided data

...
# Arguments
- `data::Tuple{Vector{BigFloat}, Vector{BigFloat}}`: Pressure dataset and geodetic altitude samples to regress against.

...
"""
function regress_poly7(data::Tuple{Vector{BigFloat}, Vector{BigFloat}})
    α_0::Vector{BigFloat} = BigFloat.(zeros(7))
    my_problem = NonlinearLeastSquaresProblem(NonlinearFunction(poly7_res), α_0, data)
    my_solution = solve(my_problem)
    return my_solution.u
end


"""
    tsos_poly7(x::Vector{T}) where T <: Real

Find the sum of squares of the given dataset

...
# Arguments
- `x::Vector{T}`: Dataset to find the total sum-of-squares of

...
"""
function tsos_poly7(x::Vector{T}) where T <: Real
    x̄ = mean(x)
    buffer = x .- x̄
    return dot(buffer, buffer)
end


"""
    rsq_poly7(data::Tuple{Vector{T}, Vector{T}}) where T<: Real

Find the coefficient of determination of the given dataset

...
# Arguments
- `α::Vector{Real}`: Seven-element vector containing the polynomial fit coefficients from 1st to 7th order
- `data::Tuple{Vector{BigFloat}, Vector{BigFloat}}`: Pressure dataset and geodetic altitude samples to regress against.

...
"""
function rsq_poly7(data::Tuple{Vector{T}, Vector{T}}) where T<: Real
    t_true, y_true = data
    res = poly7_res(Float32.(α), data)
    return 1.0 - dot(res, res) * tsos_poly7(y_true)^-1
end


function logatm7(h)
    local logp
    if h < h_base
        logp = m_base * (h - h_base) + b_base
    elseif h > h_top
        logp = m_top * (h - h_top) + b_top
    else
        logp = poly7(h, α)
    end
    return logp
end
function logatm7(h::Float32)
    local logp::Float32
    if h < Float32(h_base)
        logp = Float32(m_base) * (h - Float32(h_base)) + Float32(b_base)
    elseif h > Float32(h_top)
        logp = Float32(m_top) * (h - h_top) + Float32(b_top)
    else
        logp = poly7(h, Float32.(α))
    end
    return logp
end


setprecision(128)
SpaceIndices.init()


h_base::BigFloat = -50.0
h_top::BigFloat  = 50.0 * 10.0^3
h::Vector{BigFloat}   = collect(range(h_base, h_top, length = 2^9))
p::Vector{BigFloat}   = true_log_pressure.(h) .- true_log_pressure(0.0)
p′::Vector{BigFloat}  = true_log_deriv_pressure.(h)
p′′::Vector{BigFloat} = true_log_2nd_deriv_pressure.(h)
α::Vector{BigFloat}   = regress_poly7((h, p))
m_base::BigFloat      = true_log_deriv_pressure(h_base)
m_top::BigFloat       = true_log_deriv_pressure(h_top)
b_base::BigFloat      = true_log_pressure(h_base) .- true_log_pressure(0.0)
b_top::BigFloat       = true_log_pressure(h_top) .- true_log_pressure(0.0)
P_0::BigFloat         = pressure_nrlmsise00(0.0)


h_plot::Vector{BigFloat}   = collect(range(-50.0, 50.0 * 10.0^3, length = 2^9))
p_plot::Vector{BigFloat}   = true_log_pressure.(h_plot) .- true_log_pressure(0.0)
p_interp::Vector{BigFloat} = BigFloat.(logatm7.(Float32.(h_plot)))
rsq::BigFloat = rsq_poly7((h_plot, p_plot))
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
save(normpath(joinpath((@__FILE__), raw"..\filterbarometerfunctiondiagnostic.png")), fig, size = (3000, 2000))
fig
=#


save(normpath(joinpath((@__FILE__), raw"..\FilterBarometerData.jld2")), Dict(
    "h_base" => h_base,
    "h_top"  => h_top,
    "m_base" => m_base,
    "m_top"  => m_top,
    "b_base" => b_base,
    "b_top"  => b_top,
    "α"      => α,
    "P_0"    => P_0,
    "rsq"    => rsq,
    "rms_rel_error"   => ε_lin_rms,
    "worst_rel_error" => ε_lin_max
))


save(normpath(joinpath((@__FILE__), raw"..\FilterBarometerDataF32.jld2")), Dict(
    "h_base" => Float32(h_base),
    "h_top"  => Float32(h_top),
    "m_base" => Float32(m_base),
    "m_top"  => Float32(m_top),
    "b_base" => Float32(b_base),
    "b_top"  => Float32(b_top),
    "α"      => Float32.(α),
    "P_0"    => Float32(P_0)
))


#=
p′::Vector{BigFloat}  = true_log_deriv_pressure.(h)
p′′::Vector{BigFloat} = true_log_2nd_deriv_pressure.(h)
G::Vector{BigFloat}   = p′ ./ p
H::Vector{BigFloat}   = G.^2 - p′′ ./ p
n::Int = 5
ref_sol = (n - 1) * (n * H - G.^2)
h_f32::Vector{Float32}     = Float32.(h)
p_f32::Vector{Float32}     = Float32.(p)
p′_f32::Vector{Float32}    = Float32.(p′)
p′′_f32::Vector{Float32}   = Float32.(p′′)
G_f32::Vector{Float32}     = p′_f32 ./ p_f32
H_f32::Vector{Float32}     = G_f32.^2 - p′′_f32 ./ p_f32
f32_sol_1::Vector{Float32} = (n - 1) * (n * H_f32 - G_f32.^2)
f32_sol_2::Vector{Float32} = (n - 1) * ((n - 1) * G_f32.^2 - n * (p′′_f32 ./ p_f32))
ε_1::Vector{BigFloat} = (BigFloat.(sqrt.(f32_sol_1)) - sqrt.(ref_sol)) ./ sqrt.(ref_sol)
ε_2::Vector{BigFloat} = (BigFloat.(sqrt.(f32_sol_2)) - sqrt.(ref_sol)) ./ sqrt.(ref_sol)
# numerical errors incurred by calculating the denominator of a in single precision floating point arithmetic are at worst on the order of 10^-6 and on average on the order of 10^-9
# assuming the derivative and second derivative are correct
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, h, ε_1, label = "method 1")
lines!(ax, h, ε_2, label = "method 2")
Legend(fig, ax)
fig
# cannot evaluate a in Laguerre's method under approx. -2800 m geodetic altitude because G^2 grows larger than nH below this point
# this puts in the attractor for a complex root, not the strictly real root we seek
=#

#=
fig = Figure()
ax = Axis(fig[1, 1], title = "Test Case", xlabel = "Altitude (metres)", ylabel = "Pressure (Pascal)")
plot_samples::Vector{BigFloat} = collect(range(-10.0^3, 1.5 * 10.0^4, length = 2^8))
actual_values::Vector{BigFloat} = get_pressure.(35.34788, -117.80683, plot_samples)
lines!(ax, plot_samples, actual_values, label = "Raw Data")
lines!(ax, plot_samples, lowalt_pressure.(plot_samples, (low_coeffs, ), (h_base_low, ), (h_top_low, )), label = "Fit")
Legend(fig, ax)
fig
=#

