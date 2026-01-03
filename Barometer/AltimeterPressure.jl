begin
# interval in geodetic altitude over which the interpolant is valid
    local h_base::Float32 = -1000.0f0
    local h_ceil::Float32 = 12000.0f0
    # slope of log normalized pressure at either extrema of the valid interval
    local m_base::Float32 = -0.00011841112f0
    local m_ceil::Float32 = -0.00015328368f0
    # value of the log normalized pressure at either extrema of the valid interval
    local b_base::Float32 = 0.11881527f0
    local b_ceil::Float32 = -1.6251616f0
    # leading coefficient of the derivative of the atmospheric pressure with respect to geodetic altitude at either extrema of the valid interval
    local C_base::Float32 = -11.912918f0
    local C_ceil::Float32 = -15.421322f0
    # atmospheric pressure at the surface of the reference ellipsoid
    local P_0::Float32 = 100606.414f0
    # fifth order polynomial interpolant coefficients
    local β::Vector{Float32} = [
        -0.00011933408;
        -6.295912f-10;
        -1.06790716f-13;
        3.986928f-18;
        -2.5322159f-24
    ]


    """
        altimeter_logP(h::Float32)

    Calculates the logarithm of the normalized ambient atmospheric pressure as a function of the geodetic altitude over the WGS 84 reference ellipsoid. 5th-order polynomial interpolant used. 

    If the provided geodetic altitude exceeds the bounds in which the polynomial interpolant is valid, backup linear fits are used instead.

    ...
    # Arguments
    - `h::Float32`: Geodetic altitude above the WGS 84 reference ellipsoid (metres)

    ...
    """
    function altimeter_logP(h::Float32)
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
            logP = (β[1] * h) + (β[2] * h2) + (β[3] * h3) + (β[4] * h4) + (β[5] * h5)
        end
        return logP
    end


    """
        altimeter_P(h::Float32)

    Calculates the ambient atmospheric pressure (in Pascals) as a function of the geodetic altitude over the WGS 84 reference ellipsoid.

    The logarithm of the pressure is interpolated by a 5th-order polynomial. If the provided geodetic altitude exceeds the bounds in which the polynomial interpolant is valid, backup linear fits are used instead.

    ...
    # Arguments
    - `h::Float32`: Geodetic altitude above the WGS 84 reference ellipsoid (metres)

    ...
    """
    function altimeter_P(h::Float32)
        return P_0 * exp(altimeter_logP(h))
    end


    """
        altimeter_∂P∂h(h::Float32)

    Calculates the partial derivative of the static pressure with respect to the geodetic altitude (in Pascals per metre) over the WGS 84 reference ellipsoid.

    The logarithm of the pressure is interpolated by a 5th-order polynomial. If the provided geodetic altitude exceeds the bounds in which the polynomial interpolant is valid, backup linear fits are used instead.

    ...
    # Arguments
    - `h::Float32`: Geodetic altitude above the WGS 84 reference ellipsoid (metres)

    ...
    """
    function altimeter_∂P∂h(h::Float32)
        local ∂P∂h::Float32 = exp(altimeter_logP(h))
        if h < h_base
            ∂P∂h *= C_base 
        elseif h > h_ceil
            ∂P∂h *= C_ceil
        else
            local h2::Float32 = h^2
            local h3::Float32 = h * h2
            local h4::Float32 = h * h3
            ∂P∂h *= P_0 * (β[1] + (2.0 * β[2] * h) + (3.0 * β[3] * h2) + (4.0 * β[4] * h3) + (5.0 * β[5] * h4))
        end
        return ∂P∂h
    end
end


#=
using JLD2
data = load(normpath(joinpath((@__FILE__), raw"..\AltimeterBarometerDataF32.jld2")))
println(data)


using GLMakie
using ForwardDiff
h_plot = collect(range(-2000.0, 14.0 * 10.0^3, length = 2^11))
fig = Figure(size = (3000, 2000))
ax = Axis(fig[1, 1])
lines!(ax, h_plot, altimeter_∂P∂h.(Float32.(h_plot)), label = "interpolated pressure derivative")
lines!(ax, h_plot, ForwardDiff.derivative.((altimeter_P, ), h_plot), label = "autodiff derivative")
axislegend(ax, position = :lb)
fig
=#

