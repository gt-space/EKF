p_0::Float32 = 100606.414


α::Vector{Float32} = Float32.([
    # fit coefficients:
    -0.00011927925;
    -1.8147103f-10;
    -2.445637f-13;
    1.7510401f-17;
    -5.290156f-22;
    7.715311f-27;
    -4.4337637f-32;
    # backup linear fit for low altitudes:
    -0.00011927191;
    # backup linear fit for high altitudes:
    -0.00012193789;
])
p_filter::Vector{Float32} = Float32.([
    # polynomial fit region lower altitude bound nominal pressure:
    101208.18;
    # polynomial fit region upper altitude bound nominal pressure:
    96.20285
])
h_filter::Vector{Float32} = Float32.([
    # polynomial fit region lower altitude bound:
    -50.0;
    # polynomial fit region upper altitude bound:
    50.0 * 10.0^3
])
# used for derivative calculation:
dpdh_coeff_filter::Vector{Float32} = Float32.([
    -12.071292f0;
    -0.011730773f0;
])


β::Vector{Float32} = Float32.([
    # fit coefficients:
    -0.00011933408;
    -6.295912f-10;
    -1.06790716f-13;
    3.986928f-18;
    -2.5322159f-24;
    # backup linear fit for low altitudes:
    -0.00011841112;
    # backup linear fit for high altitudes:
    -0.00015328368;
])
p_lowalt = Float32.([
    # polynomial fit region lower altitude bound nominal pressure:
    113299.11;
    # polynomial fit region upper altitude bound nominal pressure:
    19807.377
])
h_lowalt = Float32.([
    # polynomial fit region lower altitude bound:
    -1.0 * 10.0^3;
    # polynomial fit region upper altitude bound:
    12.0 * 10.0^3;
])
# used for derivative calculation:
dpdh_coeff_lowalt::Vector{Float32} = Float32.([
    -13.415874f0;
    -3.0361476f0;
])


"""
    poly7(x::Float32, α::Vector{Float32})

Calculate the value of a seventh-order polynomial with no constant.

...
# Arguments
- `x::Float32`: Polynomial argument
- `α::Vector{Float32}`: Polynomial coefficients given as a seven-element vector

...
"""
function poly7(x::Float32, α::Vector{Float32})
    return α[1] * x + α[2] * x^2 + α[3] * x^3 + α[4] * x^4 + α[5] * x^5 + α[6] * x^6 + α[7] * x^7
end


"""
    poly7_derivative(x::Float32, α::Vector{Float32})

Calculate the value of the derivative of a seventh-order polynomial with no constant.

...
# Arguments
- `x::Float32`: Polynomial argument
- `α::Vector{Float32}`: Polynomial coefficients given as a seven-element vector

...
"""
function poly7_derivative(x::Float32, α::Vector{Float32})
    return α[1] + 2.0 * α[2] * x + 3.0 * α[3] * x^2 + 4.0 * α[4] * x^3 + 5.0 * α[5] * x^4 + 6.0 * α[6] * x^5 + 7.0 * α[7] * x^6
end


"""
    poly5(x::Float32, β::Vector{Float32})

Calculate the value of a fifth-order polynomial with no constant.

...
# Arguments
- `x::Float32`: Polynomial argument
- `β::Vector{Float32}`: Polynomial coefficients given as a five-element vector

...
"""
function poly5(x::Float32, β::Vector{Float32})
    return β[1] * x + β[2] * x^2 + β[3] * x^3 + β[4] * x^4 + β[5] * x^5
end


"""
    poly5_derivative(x::Float32, α::Vector{Float32})

Calculate the value of the derivative of a fifth-order polynomial with no constant.

...
# Arguments
- `x::Float32`: Polynomial argument
- `β::Vector{Float32}`: Polynomial coefficients given as a five-element vector

...
"""
function poly5_derivative(x::Float32, β::Vector{Float32})
    return β[1] + 2.0 * β[2] * x + 3.0 * β[3] * x^2 + 4.0 * β[4] * x^3 + 5.0 * β[5] * x^4
end


"""
    filter_pressure(h::Float32)

Calculates the approximate static pressure above a location (in Pascals), assuming `log(pressure)` follows a seventh-order interpolating polynomial in height defined by user-provided coefficients.

If the provided altitude exceeds the bounds where the polynomial interpolant is valid, backup linear fits for `log(pressure) `are used instead.

...
# Arguments
- `h::Float32`: Geodetic altitude above WGS 84 reference ellipsoid

...
"""
function filter_pressure(h::Float32)
    # in C just define α[end] & α[end - 1] as a standalone constants rather than vector elements
    # same with p_filter & h_filter
    local p::Float32
    if h < h_filter[1]
        p = p_filter[1] * exp(α[end - 1] * (h - h_filter[1]))
    elseif h > h_filter[2]
        p = p_filter[2] * exp(α[end] * (h - h_filter[2]))
    else
        p = p_0 * exp(poly7(h, α[1:7]))
    end
    return p
end


"""
    filter_pressure(h::Float32)

Calculates the approximate partial derivative of the static pressure with respect to WGS 84 geodetic altitude above a location (in Pascals per metre), assuming `log(pressure)` follows a seventh-order interpolating polynomial in height defined by user-provided coefficients.

If the provided altitude exceeds the bounds where the polynomial interpolant is valid, backup linear fits for `log(pressure) `are used instead.

...
# Arguments
- `h::Float32`: Geodetic altitude above WGS 84 reference ellipsoid

...
"""
function filter_∂p∂h(h::Float32)
    local ∂p∂h::Float32
    if h < h_filter[1]
        ∂p∂h = dpdh_coeff_filter[1] * exp(α[end - 1] * (h - h_filter[1]))
    elseif h > h_filter[2]
        ∂p∂h = dpdh_coeff_filter[2] * exp(α[end] * (h - h_filter[2]))
    else
        ∂p∂h = (p_0 * poly7_derivative(h, α[1:7])) * exp(poly7(h, α[1:7]))
    end
    return ∂p∂h
end


"""
    lowalt_pressure(h::Float32)

Calculates the approximate static pressure in Pascals above a location, assuming `log(pressure)` follows a fifth-order interpolating polynomial in height defined by user-provided coefficients.

If the provided altitude exceeds the bounds where the polynomial interpolant is valid, backup linear fits for `log(pressure) `are used instead.

Use for low altitudes.

...
# Arguments
- `h::Float32`: Geodetic altitude above WGS 84 reference ellipsoid

...
"""
function lowalt_pressure(h::Float32)
    # in C just define β[end] & β[end - 1] as a standalone constants rather than vector elements
    # same with p_lowalt & h_lowalt
    local p::Float32
    if h < h_lowalt[1]
        p = p_lowalt[1] * exp(β[end - 1] * (h - h_lowalt[1])) 
    elseif h > h_lowalt[2]
        p = p_lowalt[2] * exp(β[end] * (h - h_lowalt[2]))
    else
        p = p_0 * exp(poly5(h, β[1:5]))
    end
    return p
end


"""
    lowalt_∂p∂h(h::Float32)

Calculates the approximate partial derivative of the static pressure with respect to WGS 84 geodetic altitude above a location (in Pascals per metre), assuming `log(pressure)` follows a fifth-order interpolating polynomial in height defined by user-provided coefficients.

If the provided altitude exceeds the bounds where the polynomial interpolant is valid, backup linear fits for `log(pressure) `are used instead.

...
# Arguments
- `h::Float32`: Geodetic altitude above WGS 84 reference ellipsoid

...
"""
function lowalt_∂p∂h(h::Float32)
    local ∂p∂h::Float32
    if h < h_lowalt[1]
        ∂p∂h = dpdh_coeff_lowalt[1] * exp(β[end - 1] * (h - h_lowalt[1]))
    elseif h > h_lowalt[2]
        ∂p∂h = dpdh_coeff_lowalt[2] * exp(β[end] * (h - h_lowalt[2]))
    else
        ∂p∂h = (p_0 * poly5_derivative(h, β[1:5])) * exp(poly5(h, β[1:5]))
    end
    return ∂p∂h
end


#=
using LinearAlgebra
using GLMakie
using JLD2
using SatelliteToolboxAtmosphericModels
using ForwardDiff


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


function foo(h)
    return get_pressure(35.34788, -117.80683, h)
end


end


SpaceIndices.init()


h_samples = Float32.(collect(range(-2.0 * 10.0^3, 20.0 * 10.0^3, length = 2^13)))
h_samples_2 = Float32.(collect(range(-2.0 * 10.0^3, 20.0 * 10.0^3, length = 2^9)))


p_truth = ForwardDiff.derivative.((foo, ), collect(range(-2.0 * 10.0^3, 20.0 * 10.0^3, length = 2^9)))
p_fit = lowalt_∂p∂h.(h_samples)


fig = Figure()
ax = Axis(fig[1, 1], title = "Altimeter-Mode Pressure Function Derivative", ylabel = "pressure derivative (Pascal per metre)", xlabel = "geodetic altitude above WGS 84 reference ellipsoid (metres)")
lines!(ax, h_samples, p_fit, label = "Polynomial Fit")
scatter!(ax, h_samples_2, p_truth, markersize = 2.0, color = :red, label = "Model Samples")
axislegend(ax, position = :rb)
fig
save(normpath(joinpath((@__FILE__), raw"..\altimeter_pressure_derivative_test_plot.png")), fig, resolution = (2000, 1000), px_per_unit = 5)
fig


h_samples = Float32.(collect(range(-10.0^3, 100.0 * 10.0^3, length = 2^13)))
h_samples_2 = Float32.(collect(range(-10.0^3, 100.0 * 10.0^3, length = 2^9)))


p_truth = ForwardDiff.derivative.((foo, ), collect(range(-10.0^3, 100.0 * 10.0^3, length = 2^9)))
p_fit = filter_∂p∂h.(h_samples)


fig = Figure()
ax = Axis(fig[1, 1], title = "Filter-Mode Pressure Function Derivative", ylabel = "pressure derivative (Pascal per metre)", xlabel = "geodetic altitude above WGS 84 reference ellipsoid (metres)")
lines!(ax, h_samples, p_fit, label = "Polynomial Fit")
scatter!(ax, h_samples_2, p_truth, markersize = 2.0, color = :red, label = "Model Samples")
axislegend(ax, position = :rb)
fig
save(normpath(joinpath((@__FILE__), raw"..\pressure_derivative_test_plot.png")), fig, resolution = (2000, 1000), px_per_unit = 5)
=#


#=
using LinearAlgebra
using GLMakie
using JLD2
using SatelliteToolboxAtmosphericModels
include("Environment.jl")
SpaceIndices.init()


h_samples = Float32.(collect(range(-2.0 * 10.0^3, 20.0 * 10.0^3, length = 2^13)))
h_samples_2 = Float32.(collect(range(-2.0 * 10.0^3, 20.0 * 10.0^3, length = 2^9)))


function reference_pressure(h::Float32)
    return get_pressure(35.34788, -117.80683, h)
end


p_truth = log.(reference_pressure.(h_samples_2))
p_fit = log.(lowalt_pressure.(h_samples))


fig = Figure()
ax = Axis(fig[1, 1], title = "Altimeter-Mode Pressure Function", xlabel = "geodetic altitude above WGS 84 reference ellipsoid (metres)", ylabel = "natural log pressure (Pascal)")
lines!(ax, h_samples, p_fit, label = "Polynomial Fit")
scatter!(ax, h_samples_2, p_truth, markersize = 2.0, color = :red, label = "Model Samples")
axislegend(ax, position = :rt)
fig
save(normpath(joinpath((@__FILE__), raw"..\altimeter_test_log_plot.png")), fig, resolution = (2000, 1000), px_per_unit = 5)
fig


h_samples = Float32.(collect(range(-10.0^3, 100.0 * 10.0^3, length = 2^13)))
h_samples_2 = Float32.(collect(range(-10.0^3, 100.0 * 10.0^3, length = 2^9)))


function reference_pressure(h::Float32)
    return get_pressure(35.34788, -117.80683, h)
end


p_truth = log.(reference_pressure.(h_samples_2))
p_fit = log.(filter_pressure.(h_samples))


fig = Figure()
ax = Axis(fig[1, 1], title = "Filter-Mode Pressure Function", ylabel = "natural log pressure (Pascal)", xlabel = "geodetic altitude above WGS 84 reference ellipsoid (metres)")
lines!(ax, h_samples, p_fit, label = "Polynomial Fit")
scatter!(ax, h_samples_2, p_truth, markersize = 2.0, color = :red, label = "Model Samples")
axislegend(ax, position = :rt)
fig
save(normpath(joinpath((@__FILE__), raw"..\pressure_function_test_log_plot2.png")), fig, resolution = (2000, 1000), px_per_unit = 5)
fig
=#


















