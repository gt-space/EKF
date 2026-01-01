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
p_filter = Float32.([
    # polynomial fit region lower altitude bound nominal pressure:
    101208.18;
    # polynomial fit region upper altitude bound nominal pressure:
    96.20285
])
h_filter = Float32.([
    # polynomial fit region lower altitude bound:
    -50.0;
    # polynomial fit region upper altitude bound:
    50.0 * 10.0^3
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
    filter_pressure(h::Float32)

Calculates the approximate static pressure in Pascals above a location, assuming `log(pressure)` follows a seventh-order interpolating polynomial in height defined by user-provided coefficients.

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


using LinearAlgebra
using GLMakie
using JLD2
using SatelliteToolboxAtmosphericModels
include("Environment.jl")
SpaceIndices.init()


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
axislegend(ax, position = :lb)
fig
save(normpath(joinpath((@__FILE__), raw"..\pressure_function_test_log_plot.png")), fig, resolution = (2000, 1000), px_per_unit = 5)
fig
=#


















