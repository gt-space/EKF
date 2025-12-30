begin


### Constants
local tkjyr::Float64 = 10.0^4 * 365.25 * 86400.0 # ten thousand Julian years to seconds
local ω_Earth::Float64 = 7.292115 * 10.0^-5      # Earth angular velocity
local r_syn::Float64 = 1.00273781191135448       # ratio of synodic to sidereal period of Earth
local T_syn::Float64 = 86400.002                 # synodic period of Earth
# WGS84 reference ellipsoid properties
local a_Earth::Float64 = 6378137.0                            # semi-major axis
local f_Earth::Float64 = (298.257223563)^-1                   # flattening
local b_Earth::Float64 = a_Earth * (1.0 - f_Earth)            # semi-minor axis
local e_first::Float64 = sqrt(1.0 - (b_Earth * a_Earth^-1)^2) # first eccentricity
# secular terms from Laskar's NGT (2011)
local ψ_sec::Vector{Float64} = [0.0; 502909.66 * π * 648000.0^-1; 11119.71 * π * 648000.0^-1; 77.32 * π * 648000.0^-1; -2353.16 * π * 648000.0^-1; -180.55 * π * 648000.0^-1; 174.51 * π * 648000.0^-1; 130.95 * π * 648000.0^-1; 24.24 * π * 648000.0^-1; -47.59 * π * 648000.0^-1; -8.66 * π * 648000.0^-1]
local Θ_sec::Vector{Float64} = [0.0; -4680.93 * π * 648000.0^-1; -1.55 * π * 648000.0^-1; 1999.25 * π * 648000.0^-1; -51.38 * π * 648000.0^-1; -249.67 * π * 648000.0^-1; -39.05 * π * 648000.0^-1; 7.12 * π * 648000.0^-1; 27.87 * π * 648000.0^-1; 5.79 * π * 648000.0^-1; 2.45 * π * 648000.0^-1]
# spectral terms from Dehant & Matthews (2018) Precession, Nutation and Wobble of the Earth
local ψ_ip::Vector{Float64}  = [-17208.0 * 10.0^-3 * π * 648000.0^-1; 207.0 * 10.0^-3 * π * 648000.0^-1; 148.0 * 10.0^-3 * π * 648000.0^-1; -1317.0 * 10.0^-3 * π * 648000.0^-1; -228.0 * 10.0^-3 * π * 648000.0^-1]
local ψ_oop::Vector{Float64} = [3.3 * 10.0^-3 * π * 648000.0^-1; 0.1 * 10.0^-3 * π * 648000.0^-1; 1.2 * 10.0^-3 * π * 648000.0^-1; -1.4 * 10.0^-3 * π * 648000.0^-1; 0.3 * 10.0^-3 * π * 648000.0^-1]
local Θ_ip::Vector{Float64}  = [9205.0 * 10.0^-3 * π * 648000.0^-1; -90.0 * 10.0^-3 * π * 648000.0^-1; 7.0 * 10.0^-3 * π * 648000.0^-1; 573.0 * 10.0^-3 * π * 648000.0^-1; 98.0 * 10.0^-3 * π * 648000.0^-1]
local Θ_oop::Vector{Float64} = [1.5 * 10.0^-3 * π * 648000.0^-1; 0.0 * 10.0^-3 * π * 648000.0^-1; -0.2 * 10.0^-3 * π * 648000.0^-1; -0.5 * 10.0^-3 * π * 648000.0^-1; 0.1 * 10.0^-3 * π * 648000.0^-1]
local ω_nut::Vector{Float64} = [(-18.6 * 365.25 * 86400.0)^(-1); (-9.3 * 365.25 * 86400.0)^(-1); (365.25 * 86400.0)^(-1); (0.5 * 365.25 * 86400.0)^(-1); (13.66 * 86400.0)^(-10)]
local χ_nut::Vector{Float64} = [0.0; 0.0; 0.0; 0.0; 0.0]


"""
    lla2ecef(φ::Float64, λ::Float64, h::Float64)

Returns the ECEF reference frame position vector corresponding to the given geodetic LLA state.

...
# Arguments
- `φ::Float64`: Geodetic latitude
- `λ::Float64`: Geodetic longitude
- `h::Float64`: Geodetic altitude

...
"""
function lla2ecef(φ::Float64, λ::Float64, h::Float64)
    local N::Float64 = a_Earth * √(1.0 - (e_first^2 * sin(φ)^2))^-1
    return r_ECEF::Vector{Float64} = [
    (N + h) * cos(φ) * cos(λ);
    (N + h) * cos(φ) * sin(λ);
    ((b_Earth^2 * a_Earth^-2) * N + h) * sin(φ)
    ]
end


"""
    ecef2lla(r::Vector{Float64}, n_max::Integer = 128, ε::Float64 = 10.0^-13)

Returns the geodetic LLA state (geodetic latitude, geodetic longitude, and geodetic altitude) corresponding to the given ECEF reference frame position vector.

...
# Arguments
- `r::Vector{Float64}`: ECEF reference frame position vector.
- `n_max::Integer`: Maximum allowed number of Newton-Raphson solves. `128` by default.
- `ε::Float64`: Allowed truncation error. `10.0^-13` by default.

...
"""
function ecef2lla(r::Vector{Float64}, n_max::Integer = 128, ε::Float64 = 10.0^-13)
    local p::Float64 = √(r[1]^2 + r[2]^2)
    local z_norm::Float64 = abs(r[3])
    local φ::Float64
    local λ::Float64
    local h::Float64
    if p == 0.0
        if r[3] == 0.0
            φ = 0.0
            λ = 0.0
            h = -b_Earth
        else
            if r[3] > 0
                φ = 0.5π
            else
                φ = -0.5π
            end
            λ = 0.0
            h = z_norm - b_Earth
        end
    elseif r[3] == 0.0
        φ = 0.0
        λ = atan(r[2], r[1])
        h = p - a_Earth
    else
        # in the non-trivial case, we are essentially solving Bowring's formula, but parameterized into a more efficient quartic equation
        # this saves computation time by avoiding unnecessary expensive calls to transcendental functions
        # to avoid the issue of the Newton-Raphson solver converging on the incorrect root, we must first initialize the iteration variable, t, for a number of different circumstances: 
        local c::Float64 = a_Earth * e_first^2
        local e′::Float64 = √(1.0 - e_first^2)
        local z′::Float64 = e′ * z_norm
        local u::Float64 = 2.0 * (z′ - c)
        local v::Float64 = 2.0 * (z′ + c)
        local t_M::Float64 = (c - z′) * p^-1
        local f_M::Float64 = p * t_M^4 + u * t_M^3 + v * t_M - p
        local t_1::Float64 = (p - c + z′) * (p - c + 2.0 * z′)^-1
        local t_0::Float64 = p * (z′ + c)^-1
        local t::Float64
        local Δt::Float64
        if t_M <= 0
            t = t_1
        elseif t_M >= 1 
            t = t_0
        elseif f_M >= 0
            t = t_0
        else
            t = t_1
        end
        for i = 1:n_max
            Δt = (p - (p * t^4 + u * t^3 + v * t)) * (4.0 * p * t^3 + 3.0 * u * t^2 + v)^-1
            t += Δt
            if abs(Δt) <= ε
                break
            end
        end
        if r[3] > 0
        φ = atan(1.0 - t^2, 2.0 * e′ * t)
        else
            φ = -atan(1.0 - t^2, 2.0 * e′ * t)
        end
        h = ((2.0 * p * e′ * t) + (z_norm * (1.0 - t^2)) - (a_Earth * e′ * (1.0 + t^2))) * √((1.0 + t^2)^2 - (4.0 * e_first^2 * t^2))^-1
        λ = atan(r[2], r[1])
    end
    if λ < 0
        λ = 2.0π + λ
    end
    return φ, λ, h
end


"""
    ecef2enu(φ::Float64, λ::Float64)

Generates the transformation matrix from the ECEF reference frame to the ENU reference frame.

...
# Arguments
- `φ::Float64`: Geodetic latitude of the ENU reference frame origin
- `λ::Float64`: Geodetic longitude of the ENU reference frame origin

...
"""
function ecef2enu(φ::Float64, λ::Float64)
    return ecef2enu_mat::Matrix{Float64} = [
        -sin(λ)        cos(λ)       0.0;
        -sin(φ)cos(λ) -sin(φ)sin(λ) cos(φ);
        cos(φ)cos(λ)   cos(φ)sin(λ) sin(φ)
    ]
end


"""
    enu2ecef(φ::Float64, λ::Float64)

Generates the transformation matrix from the ENU reference frame to the ECEF reference frame.

...
# Arguments
- `φ::Float64`: Geodetic latitude of the ENU reference frame origin
- `λ::Float64`: Geodetic longitude of the ENU reference frame origin

...
"""
function enu2ecef(φ::Float64, λ::Float64)
    return enu2ecef_mat::Matrix{Float64} = [
        -sin(λ) -sin(φ)cos(λ) cos(φ)cos(λ);
        cos(λ)  -sin(φ)sin(λ) cos(φ)sin(λ);
        0.0      cos(φ)       sin(φ)
    ]
end


"""
    get_ωα_ECEF(t::Float64, Δt_iJ2k::Float64, Δt_iY::Float64)

Returns the angular velocity and angular acceleration of the ECEF frame relative to the ECI frame for a given time.

...
# Arguments
- `t::Float64`: Time elapsed from reference time in seconds.
- `Δt_iJ2k::Float64`: Time elapsed between reference time and the J2000 epoch.
- `Δt_iY::Float64`: Time elapsed between reference time and the start of the year.

...
"""
function get_ωα_ECEF(t::Float64, Δt_iJ2k::Float64, Δt_iY::Float64)
    # reduced time interval between current time and J2000 epoch
    local t̄::Float64 = (t + Δt_iJ2k) * tkjyr^-1
    # secular precession
    local secpol::Vector{Float64}   = [1.0; t̄; t̄^2; t̄^3; t̄^4; t̄^5; t̄^6; t̄^7; t̄^8; t̄^9; t̄^10]
    local secpol′::Vector{Float64}  = [0.0; 1.0; 2.0t̄; 3.0t̄^2; 4.0t̄^3; 5.0t̄^4; 6.0t̄^5; 7.0t̄^6; 8.0t̄^7; 9.0t̄^8; 10.0t̄^9]
    local secpol′′::Vector{Float64} = [0.0; 0.0; 2.0; 6.0t̄; 12.0t̄^2; 20.0t̄^3; 30.0t̄^4; 42.0t̄^5; 56.0t̄^6; 72.0t̄^7; 90.0t̄^8]
    # spectral precession
    local Ξ_vec::Vector{Float64} = ((t + Δt_iJ2k) * ω_nut) + χ_nut
    local odds::Vector{Float64}  = sin.(Ξ_vec)
    local evens::Vector{Float64} = cos.(Ξ_vec)
    # combine terms with their respective coefficients to get angular offset and angular rates
    local ψ::Float64 = dot(ψ_sec, secpol) + dot(ψ_ip, odds) + dot(ψ_oop, evens)
    local θ::Float64 = dot(Θ_sec, secpol) + dot(Θ_ip, evens) + dot(Θ_oop, odds)
    local φ::Float64 = ω_Earth * (t + Δt_iY - ((T_syn * r_syn^-1) * floor(t + Δt_iY * r_syn * T_syn^-1)))
    local ψ̇::Float64 = dot(ψ_sec, secpol′) * tkjyr^-1 + dot(ω_nut .* ψ_ip, evens) - dot(ω_nut .* ψ_oop, odds)
    local θ̇::Float64 = dot(Θ_sec, secpol′) * tkjyr^-1 - dot(ω_nut .* Θ_ip, odds) + dot(ω_nut .* Θ_oop, evens)
    local ψ̇̇::Float64 = dot(ψ_sec, secpol′′) * tkjyr^-2 - dot(ω_nut.^2 .* ψ_ip, odds) - dot(ω_nut.^2 .* ψ_oop, evens)
    local θ̇̇::Float64 = dot(Θ_sec, secpol′′) * tkjyr^-2 - dot(ω_nut.^2 .* Θ_ip, evens) - dot(ω_nut.^2 .* Θ_oop, odds)
    # compute vectors
    local cθ::Float64 = cos(θ)
    local sθ::Float64 = sin(θ)
    local cψ::Float64 = cos(ψ)
    local sψ::Float64 = sin(ψ)
    return ω::Vector{Float64} = [
        (ω_Earth * sθ * sψ) + (θ̇ * cψ);
        (-ω_Earth * sθ * cψ) + (θ̇ * sψ);
        (ω_Earth * cθ) + ψ̇ 
    ],
    α::Vector{Float64} = [
        (θ̇ * ω_Earth * cθ * sψ) + (θ̇̇ * cψ) + (ψ̇ * ((ω_Earth * sθ * cψ) - (θ̇ * sψ)));
        (-θ̇ * ω_Earth * cθ * cψ) + (θ̇̇ * sψ) + (ψ̇ * ((ω_Earth * sθ * sψ) + (θ̇ * cψ)));
        (-θ̇ * ω_Earth * sθ) + ψ̇̇
    ]
end


#=
using AstroTime
#year month day hour minute second
current_epoch = Epoch{UniversalTime}(2026, 1, 10, 12, 00, 0.0)
reference_epoch = UT1Epoch(0.0seconds, origin = :j2000)
year_epoch = Epoch{UniversalTime}(2000, 1, 1, 0, 0, 0.0)
Δt_iJ2k = value(seconds(current_epoch - reference_epoch))
Δt_iY = value(seconds(current_epoch - year_epoch))
get_ωα_ECEF(0.0, Δt_iJ2k, Δt_iY)
=#


end





