begin
    # interval in geodetic altitude over which the interpolant is valid
    local h_base::Float32 = -50.0f0
    local h_ceil::Float32 = 50000.0f0
    # slope of log normalized pressure at either extrema of the valid interval
    local m_base::Float32 = -0.00011927191f0
    local m_ceil::Float32 = -0.00012193789f0
    # value of the log normalized pressure at either extrema of the valid interval
    local b_base::Float32 = 0.0059635397f0
    local b_ceil::Float32 = -6.9525123f0
    # leading coefficient of the derivative of the atmospheric pressure with respect to geodetic altitude at either extrema of the valid interval
    local C_base::Float32 = -11.999518f0
    local C_ceil::Float32 = -12.267735f0
    # atmospheric pressure at the surface of the reference ellipsoid
    local P_0::Float32 = 100606.414f0
    # seventh order polynomial interpolant coefficients
    local α::Vector{Float32} = [
        -0.00011927925;
        -1.8147103f-10;
        -2.445637f-13;
        1.7510401f-17;
        -5.290156f-22;
        7.715311f-27;
        -4.4337637f-32
    ]

    
    """
        filter_logP(h::Float32)

    Calculates the logarithm of the normalized ambient atmospheric pressure as a function of the geodetic altitude over the WGS 84 reference ellipsoid. 7th-order polynomial interpolant used. 

    If the provided geodetic altitude exceeds the bounds in which the polynomial interpolant is valid, backup linear fits are used instead.

    ...
    # Arguments
    - `h::Float32`: Geodetic altitude above the WGS 84 reference ellipsoid (metres)

    ...
    """
    function filter_logP(h::Float32)
        local logP::Float32
        # you could probably optimize this further by splitting the conditional blocks into their own unique functions
        if h < h_base
            logP = m_base * (h - h_base) + b_base
        elseif h > h_ceil
            logP = m_ceil * (h - h_ceil) + b_ceil
        else
            local h2::Float32 = h^2
            local h3::Float32 = h * h2
            local h4::Float32 = h * h3
            local h5::Float32 = h * h4
            local h6::Float32 = h * h5
            local h7::Float32 = h * h6
            logP = (α[1] * h) + (α[2] * h2) + (α[3] * h3) + (α[4] * h4) + (α[5] * h5) + (α[6] * h6) + (α[7] * h7)
        end
        return logP
    end


    """
        filter_P(h::Float32)

    Calculates the ambient atmospheric pressure (in Pascals) as a function of the geodetic altitude over the WGS 84 reference ellipsoid.

    The logarithm of the pressure is interpolated by a 7th-order polynomial. If the provided geodetic altitude exceeds the bounds in which the polynomial interpolant is valid, backup linear fits are used instead.

    ...
    # Arguments
    - `h::Float32`: Geodetic altitude above the WGS 84 reference ellipsoid (metres)

    ...
    """
    function filter_P(h::Float32)
        return P_0 * exp(filter_logP(h))
    end


    """
        filter_∂P∂h(h::Float32)

    Calculates the partial derivative of the static pressure with respect to the geodetic altitude (in Pascals per metre) over the WGS 84 reference ellipsoid.

    The logarithm of the pressure is interpolated by a 7th-order polynomial. If the provided geodetic altitude exceeds the bounds in which the polynomial interpolant is valid, backup linear fits are used instead.

    ...
    # Arguments
    - `h::Float32`: Geodetic altitude above the WGS 84 reference ellipsoid (metres)

    ...
    """
    function filter_∂P∂h(h::Float32)
        local ∂P∂h::Float32 = exp(filter_logP(h))
        if h < h_base
            ∂P∂h *= C_base 
        elseif h > h_ceil
            ∂P∂h *= C_ceil
        else
            local h2::Float32 = h^2
            local h3::Float32 = h * h2
            local h4::Float32 = h * h3
            local h5::Float32 = h * h4
            local h6::Float32 = h * h5
            ∂P∂h *= P_0 * (α[1] + (2.0 * α[2] * h) + (3.0 * α[3] * h2) + (4.0 * α[4] * h3) + (5.0 * α[5] * h4) + (6.0 * α[6] * h5) + (7.0 * α[7] * h6))
        end
        return ∂P∂h
    end
end


#=
using JLD2
data = load(normpath(joinpath((@__FILE__), raw"..\FilterBarometerDataF32.jld2")))


using GLMakie
using ForwardDiff
h_plot = collect(range(-1000.0, 60.0 * 10.0^3, length = 2^11))
fig = Figure(size = (3000, 2000))
ax = Axis(fig[1, 1])
lines!(ax, h_plot, filter_∂P∂h.(Float32.(h_plot)), label = "interpolated pressure derivative")
lines!(ax, h_plot, ForwardDiff.derivative.((filter_P, ), h_plot), label = "autodiff derivative")
axislegend(ax, position = :lb)
fig
=#

