using LinearAlgebra
using AstroTime


global a::Float64 = 6378136.3   # semi-major axis
global f::Float64 = (298.257)^-1  # flattening
global b::Float64 = a * (1.0 - f) # semi-minor axis
global e::Float64 = sqrt(1.0 - (b^2 / a^2)) # first eccentricity
global ϵ::Float64 = 10.0^-12 # allowable error
global i_max = 100 # maximum allowable iterations before truncation


global ω_0::Float64 = 7.292115 * 10^-5
global r_syn::Float64 = 1.00273781191135448
global T_syn::Float64 = 86400.002
# Secular terms from Laskar's NGT
global ψ_sec::Vector{Float64} = [0.0; 502909.66 * π * 648000.0^-1; 11119.71 * π * 648000.0^-1; 77.32 * π * 648000.0^-1; -2353.16 * π * 648000.0^-1; -180.55 * π * 648000.0^-1; 174.51 * π * 648000.0^-1; 130.95 * π * 648000.0^-1; 24.24 * π * 648000.0^-1; -47.59 * π * 648000.0^-1; -8.66 * π * 648000.0^-1]
global Θ_sec::Vector{Float64} = [0.0; -4680.93 * π * 648000.0^-1; -1.55 * π * 648000.0^-1; 1999.25 * π * 648000.0^-1; -51.38 * π * 648000.0^-1; -249.67 * π * 648000.0^-1; -39.05 * π * 648000.0^-1; 7.12 * π * 648000.0^-1; 27.87 * π * 648000.0^-1; 5.79 * π * 648000.0^-1; 2.45 * π * 648000.0^-1]
# Spectral terms from Dehant & Matthews
global ψ_ip::Vector{Float64}  = [-17208.0 * 10.0^-3 * π * 648000.0^-1; 207.0 * 10.0^-3 * π * 648000.0^-1; 148.0 * 10.0^-3 * π * 648000.0^-1; -1317.0 * 10.0^-3 * π * 648000.0^-1; -228.0 * 10.0^-3 * π * 648000.0^-1]
global ψ_oop::Vector{Float64} = [3.3 * 10.0^-3 * π * 648000.0^-1; 0.1 * 10.0^-3 * π * 648000.0^-1; 1.2 * 10.0^-3 * π * 648000.0^-1; -1.4 * 10.0^-3 * π * 648000.0^-1; 0.3 * 10.0^-3 * π * 648000.0^-1]
global Θ_ip::Vector{Float64}  = [9205.0 * 10.0^-3 * π * 648000.0^-1; -90.0 * 10.0^-3 * π * 648000.0^-1; 7.0 * 10.0^-3 * π * 648000.0^-1; 573.0 * 10.0^-3 * π * 648000.0^-1; 98.0 * 10.0^-3 * π * 648000.0^-1]
global Θ_oop::Vector{Float64} = [1.5 * 10.0^-3 * π * 648000.0^-1; 0.0 * 10.0^-3 * π * 648000.0^-1; -0.2 * 10.0^-3 * π * 648000.0^-1; -0.5 * 10.0^-3 * π * 648000.0^-1; 0.1 * 10.0^-3 * π * 648000.0^-1]
global ω_nut::Vector{Float64} = [(-18.6 * 365.25 * 86400.0)^(-1); (-9.3 * 365.25 * 86400.0)^(-1); (365.25 * 86400.0)^(-1); (0.5 * 365.25 * 86400.0)^(-1); (13.66 * 86400.0)^(-10)]
global χ_nut::Vector{Float64} = [0.0; 0.0; 0.0; 0.0; 0.0]


function get_ωα_ECEF(t::Float64, Δτ_iJ2k::Float64, Δτ_iY::Float64)
    t̄::Float64 = (t + Δτ_iJ2k) * (10.0^4 * 365.25 * 86400.0)^-1
    secpol::Vector{Float64}   = [1.0; t̄; t̄^2; t̄^3; t̄^4; t̄^5; t̄^6; t̄^7; t̄^8; t̄^9; t̄^10]
    secpol′::Vector{Float64}  = [0.0; 1.0; 2.0t̄; 3.0t̄^2; 4.0t̄^3; 5.0t̄^4; 6.0t̄^5; 7.0t̄^6; 8.0t̄^7; 9.0t̄^8; 10.0t̄^9]
    secpol′′::Vector{Float64} = [0.0; 0.0; 2.0; 6.0t̄; 12.0t̄^2; 20.0t̄^3; 30.0t̄^4; 42.0t̄^5; 56.0t̄^6; 72.0t̄^7; 90.0t̄^8]
    Ξ_vec::Vector{Float64} = ((t + Δτ_iJ2k) * ω_nut) + χ_nut
    odds::Vector{Float64}  = sin.(Ξ_vec)
    evens::Vector{Float64} = cos.(Ξ_vec)


    ψ::Float64 = dot(ψ_sec, secpol) + dot(ψ_ip, odds) + dot(ψ_oop, evens)
    θ::Float64 = dot(Θ_sec, secpol) + dot(Θ_ip, evens) + dot(Θ_oop, odds)
    φ::Float64 = ω_0 * (t + Δτ_iY - ((T_syn * r_syn^-1) * floor(t + Δτ_iY * r_syn * T_syn^-1)))
    ψ̇::Float64 = dot(ψ_sec, secpol′) * (10.0^4 * 365.25 * 86400.0)^-1 + dot(ω_nut .* ψ_ip, evens) - dot(ω_nut .* ψ_oop, odds)
    θ̇::Float64 = dot(Θ_sec, secpol′) * (10.0^4 * 365.25 * 86400.0)^-1 - dot(ω_nut .* Θ_ip, odds) + dot(ω_nut .* Θ_oop, evens)
    ψ̇̇::Float64 = dot(ψ_sec, secpol′′) * (10.0^4 * 365.25 * 86400.0)^-2 - dot(ω_nut.^2 .* ψ_ip, odds) - dot(ω_nut.^2 .* ψ_oop, evens)
    θ̇̇::Float64 = dot(Θ_sec, secpol′′) * (10.0^4 * 365.25 * 86400.0)^-2 - dot(ω_nut.^2 .* Θ_ip, evens) - dot(ω_nut.^2 .* Θ_oop, odds)


    cθ::Float64 = cos(θ)
    sθ::Float64 = sin(θ)
    cψ::Float64 = cos(ψ)
    sψ::Float64 = sin(ψ)
    ω::Vector{Float64} = [
        (ω_0 * sθ * sψ) + (θ̇ * cψ);
        (-ω_0 * sθ * cψ) + (θ̇ * sψ);
        (ω_0 * cθ) + ψ̇ 
    ]
    α::Vector{Float64} = [
        (θ̇ * ω_0 * cθ * sψ) + (θ̇̇ * cψ) + (ψ̇ * ((ω_0 * sθ * cψ) - (θ̇ * sψ)));
        (-θ̇ * ω_0 * cθ * cψ) + (θ̇̇ * sψ) + (ψ̇ * ((ω_0 * sθ * sψ) + (θ̇ * cψ)));
        (-θ̇ * ω_0 * sθ) + ψ̇̇
    ]
    return ω, α
end


# geo_lla = [φ; λ; h]
function lla2ecef(φ::Float64; λ::Float64, h::Float64)
    local N::Float64 = a * √(1.0 - (e^2 * sin(φ)^2))^-1
    local r_ECEF::Vector{Float64} = [
    (N + h) * cos(φ) * cos(λ);
    (N + h) * cos(φ) * sin(λ);
    ((b^2 * a^-2) * N + h) * sin(φ)
    ]
    return r_ECEF
end


function ecef2lla(r::Vector{Float64})
    local p::Float64 = √(r[1]^2 + r[2]^2)
    local z_norm::Float64 = abs(r[3])
    local φ::Float64
    local λ::Float64
    local h::Float64
    if p == 0.0
        if r[3] == 0.0
            φ = 0.0
            λ = 0.0
            h = -b
        else
            if r[3] > 0
                φ = 0.5π
            else
                φ = -0.5π
            end
            λ = 0.0
            h = z_norm - b
        end
    elseif r[3] == 0.0
        φ = 0.0
        λ = atan(r[2], r[1])
        h = p - a
    else
        # in the non-trivial case, we are essentially solving Bowring's formula, but parameterized into a more efficient quartic equation
        # this saves computation time by avoiding unnecessary expensive calls to transcendental functions
        # to avoid the issue of the Newton-Raphson solver converging on the incorrect root, we must first initialize the iteration variable, t, for a number of different circumstances: 
        local c::Float64 = a * e^2
        local e′::Float64 = √(1.0 - e^2)
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
        for i = 1:i_max
            Δt = (p - (p * t^4 + u * t^3 + v * t)) * (4.0 * p * t^3 + 3.0 * u * t^2 + v)^-1
            t += Δt
            if abs(Δt) <= ϵ
                break
            end
        end
        if r[3] > 0
        φ = atan(1.0 - t^2, 2.0 * e′ * t)
        else
            φ = -atan(1.0 - t^2, 2.0 * e′ * t)
        end
        h = ((2.0 * p * e′ * t) + (z_norm * (1.0 - t^2)) - (a * e′ * (1.0 + t^2))) * √((1.0 + t^2)^2 - (4.0 * e^2 * t^2))^-1
        λ = atan(r[2], r[1])
    end
    if λ < 0
        λ = 2.0π + λ
    end
    geo_lla = [φ; λ; h]
end


function zyx2q(ψ::Float64, θ::Float64, φ::Float64)
    i::T = cos(ψ / 2) * cos(θ / 2) * sin(φ / 2) - sin(ψ / 2) * sin(θ / 2) * cos(φ / 2)
    j::T = cos(ψ / 2) * sin(θ / 2) * cos(φ / 2) + sin(ψ / 2) * cos(θ / 2) * sin(φ / 2)
    k::T = sin(ψ / 2) * cos(θ / 2) * cos(φ / 2) - cos(ψ / 2) * sin(θ / 2) * sin(φ / 2)
    r::T = cos(ψ / 2) * cos(θ / 2) * cos(φ / 2) + sin(ψ / 2) * sin(θ / 2) * sin(φ / 2)
    q = Quaternion(i, j, k, r)
    q /= abs(q)
    return q
end



#year month day hour minute second
current_epoch = Epoch{UniversalTime}(2026, 1, 10, 12, 00, 0.0)
reference_epoch = UT1Epoch(0.0seconds, origin = :j2000)
year_epoch = Epoch{UniversalTime}(2000, 1, 1, 0, 0, 0.0)
Δτ_iJ2k = value(seconds(current_epoch - reference_epoch))
Δτ_iY = value(seconds(current_epoch - year_epoch))
get_ωα_ECEF(0.0, Δτ_iJ2k, Δτ_iY)



