begin
    # Julia's compiler can automatically detect which local values are constant, so no local const declaration is required


    # Atmosphere constants:
    # interval in geodetic altitude over which the interpolant is valid
    local h_base::Float32 = -1000.0f0
    local h_ceil::Float32 = 12000.0f0
    # slope of log normalized pressure at either extrema of the valid interval
    local m_base::Float32 = -0.00011841112f0
    local m_ceil::Float32 = -0.00015328368f0
    # value of the log normalized pressure at either extrema of the valid interval
    local b_base::Float32 = 0.11881527f0
    local b_ceil::Float32 = -1.6251616f0
    # atmospheric pressure at the surface of the reference ellipsoid
    local P_0::Float32 = 100606.414f0
    # log atmospheric pressure at the surface of the reference ellipsoid
    local b_0::Float32 = 11.518971f0
    # fifth order polynomial interpolant coefficients
    local β::Vector{Float32} = [
        -0.00011933408;
        -6.295912f-10;
        -1.06790716f-13;
        3.986928f-18;
        -2.5322159f-24
    ]


    # Solver constants:
    # interpolant polynomial order
    local n::Int32 = 5
    # maximum allowable numerical solver error
    local ε::Float32 = 1.0f-7
    # maximum allowable number of solver iterations
    local i_max::Int32 = 2


    function Laguerre_solve(x_0::Float32, f̂::Float32)
        local x::Float32 = x_0
        local x2::Float32
        local x3::Float32
        local x4::Float32
        local x5::Float32
        local f::Float32
        local f′::Float32
        local f′′::Float32
        local G::Float32
        local H::Float32
        local λ::Float32
        local a::Float32
        for i = 1:i_max # make sure there is an interrupt
            # buffer for storing the powers in x
            x2 = x^2
            x3 = x * x2
            x4 = x * x3
            x5 = x * x4
            # objective function value:
            f = (β[1] * x) + (β[2] * x2) + (β[3] * x3) + (β[4] * x4) + (β[5] * x5) - f̂
            # this is a safety feature to avoid division by zero:
            if abs(f) < ε
                break
            end
            # objective function derivative:
            f′ = (β[1]) + (2.0 * β[2] * x) + (3.0 * β[3] * x2) + (4.0 * β[4] * x3) + (5.0 * β[5] * x4)
            # objective function second derivative:
            f′′ = (2.0 * β[2]) + (6.0 * β[3] * x) + (12.0 * β[4] * x2) + (20.0 * β[5] * x3)
            # G, H, and a are values used in Laguerre's method
            G = f′ / f
            H = G^2 - (f′′ / f)
            # λ is just a buffer to store the denominator of a
            # factoring the n over the G^2 in H in this expression is inadvisable as this approximately doubles the floating point noise and round-off errors
            λ = √((n - 1) * (n * H - G^2))
            # this upcoming one liner expression for a is identical to the following conditional block:
            # if G < 0.0
            #    a = n / (G - λ)
            # else
            #    a = n / (G + λ)
            # end
            # essentially, when calculating a, we want to minimize the step (because we are finding the closest root)
            # additionally, we want to maximize the denominator of a, which is G ± λ
            # in order to avoid catastrophic cancellation, we pick the sign of ± such that G ± λ is maximized, which of course also minimizes a
            a = G < 0.0 ? (n / (G - λ)) : (n / (G + λ))
            # update step:
            x -= a
            # if the change just made was small, we have arrived at the solution
            if abs(a) < ε
                break
            end
        end
        return x
    end


    """
        lerp(logP::Float32)

    Estimate the geodetic altitude (in metres) by lerping between the ends of the valid interval, using the logarithm of the normalized atmospheric pressure as the lerp parameter.

    ...
    # Arguments
    - `logP::Float32`: Log normalized atmospheric pressure (Pascal)

    ...
    """
    function lerp(logP::Float32)
        return (h_ceil * (b_base - logP) + h_base * (logP - b_ceil)) / (b_base - b_ceil)
    end

    
    """
        logP2alt(logP::Float32)

    Calculate the geodetic altitude (in metres) as a function of the logarithm of the normalized atmospheric pressure.

    The normalized atmospheric pressure is the pressure divided by the pressure at the surface of the WGS 84 reference ellipsoid.

    ...
    # Arguments
    - `logP::Float32`: Log normalized atmospheric pressure (Pascal)

    ...
    """
    function logP2alt(logP::Float32)
        local h::Float32
        # first two cases happen if logP is out of bounds
        # if logP is greater than the maximum interpolant pressure, we are below h_base
        if logP > b_base
            h = h_base + (logP - b_base) / m_base
        # if logP is smaller than the minimum interpolant pressure, we are above h_ceil
        elseif logP < b_ceil
            h = h_ceil + (logP - b_ceil) / m_ceil
        # otherwise we will call Laguerre's method 
        else
            # in order to initialize Laguerre's method, we need an initial guess for our altitude
            # we will use a single linsolve in order to achieve this
            h = Laguerre_solve(
                lerp(logP),
                logP
            )
        end
        return h
    end


    """
        pressure_altimeter_uncorrected(P::Float32)

    Calculate the geodetic altitude (in metres) as a function of the atmospheric pressure.

    ...
    # Arguments
    - `P::Float32`: Atmospheric pressure (Pascal)

    ...
    """
    function pressure_altimeter_uncorrected(P::Float32)
        return logP2alt(log(P) - b_0)
    end


    """
        pressure_altimeter_corrected(P::Float32, Δh::Float32)

    Calculate the geodetic altitude (in metres) as a function of the atmospheric pressure.

    ...
    # Arguments
    - `P::Float32`: Atmospheric pressure (Pascal)
    - `Δh::Float32`: Height offset used for compensation (Pascal)

    ...
    """
    function pressure_altimeter_corrected(P::Float32, Δh::Float32)
        return pressure_altimeter_uncorrected(P) - Δh
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

