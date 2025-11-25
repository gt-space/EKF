using LinearAlgebra


using JLD2


include("CoordinateTransformations.jl")


perf_dict = load(normpath(joinpath((@__FILE__), raw"..\PerformanceData.jld2")))


global μ::Float64 = 3.986004418 * 10.0^14
global J2_E::Float64 = 1.08262668 * 10.0^-3
global R_E::Float64 = 6.3710008 * 10.0^6
global g_0::Float64 = 9.80665


global CN_LUT::Vector{Float64}    = perf_dict["CN"]
global CM_LUT::Vector{Float64}    = perf_dict["CM"]
global CAon_LUT::Matrix{Float64}  = perf_dict["CAon"]
global CAoff_LUT::Matrix{Float64} = perf_dict["CAoff"]
global Δα::Float64 = π * 180.0^-1 # radians (spreadsheet has steps in 1 degree)
global ΔMach::Float64 = 0.01 # dimensionless
global Mach_samples::Vector{Float64} = collect(range(0.01, 4.5, step = 0.01))
global α_samples::Vector{Float64} = π * 180.0^-1 * collect(range(0, 15, step = 1))


global h_samples::Vector{Float64}   = perf_dict["heights"]
global ρ_LUT::Vector{Float64}       = perf_dict["densities"]
global temp_LUT::Vector{Float64}    = perf_dict["temperatures"]
global molm_LUT::Vector{Float64}    = perf_dict["molar_masses"]
global pstatic_LUT::Vector{Float64} = perf_dict["pressures"]
global γ::Float64 = 1.4 # this is only an approximation for the adiabatic index
global Δh::Float64 = 10.0^3


global m_dry_nom::Float64 = 128.48
global m_fuel_nom::Float64 = 40.25
global m_ox_nom::Float64 = 80.85
global m_gas_nom::Float64 = 15.0
global m_wet_nom::Float64 = m_dry_nom + m_fuel_nom + m_ox_nom + m_gas_nom
# simulation body frame y is positive towards nose
# z in CAD is x in simulation body frame
# y in CAD is z in simulation body frame

# y in CAD is x in simulation body frame
# z in CAD is y in simulation body frame
l_OA::Float64 = 7.23 + 0.15875
global d_AVI::Float64 = 41.03 * 0.0254
global d̂_thrust::Vector{Float64} = [0.0; 0.0; 1.0]
global r_thrust::Vector{Float64} = [0.0; 0.0; -l_OA] + [0.0; 0.0; d_AVI]
global cg_dry::Vector{Float64} = r_thrust + [0.000762; 0.00127; 4.18719]
global cg_wet::Vector{Float64} = r_thrust + [0.000508; -0.000254; 3.681]
# mdots in kg * s^-1
global mdot_fuel::Float64 = -3.2
global mdot_ox::Float64 = -1.6
# thrust in Newtons
global thrust::Float64 = 11120.554
# pressure in Pascals
global p_nozzle::Float64 = 68947.0
# areas in m^2
global S_nozzle::Float64 = 0.04047089
global S_ref::Float64 = 0.04047089


global I_dry::Matrix{Float64} = [
    3547.5412 -0.01225 -1.0125;
    -0.01225 3547.558 -2.4294;
    -1.0125 -2.4294 6.3856
] + m_dry_nom * (dot(cg_dry, cg_dry) * I(3) - cg_dry * transpose(cg_dry))
global I_wet::Matrix{Float64} = [
    4588.36 -0.0119 -0.579;
    -0.0119 4588.36 -2.194;
    -0.579 -2.194 10.1713
] + m_wet_nom * (dot(cg_wet, cg_wet) * I(3) - cg_wet * transpose(cg_wet))
global I_dry_inv = inv(I_dry)
global I_wet_inv = inv(I_wet) 


function bracket(x::Float64, samples::Vector{Float64})
    local a::Integer
    local b::Integer
    if x <= samples[1]
        a = 1
        b = 2
    elseif x >= samples[end]
        b = length(samples)
        a = b - 1
    else
        a = findlast(x .> samples)
        b = a + 1
    end
    return a, b
end


function lerp(x::Float64, Δx::Float64, x_samples::Vector{Float64}, LUT::Vector{Float64})
    a, b = bracket(x, x_samples)
    return (LUT[a] * (x_samples[b] - x) + LUT[b] * (x - x_samples[a])) * Δx^-1
end


function lerp_with_deriv(x::Float64, Δx::Float64, x_samples::Vector{Float64}, LUT::Vector{Float64})
    a, b = bracket(x, x_samples)
    return (LUT[a] * (x_samples[b] - x) + LUT[b] * (x - x_samples[a])) * Δx^-1, (x_samples[b] - x_samples[a]) * Δx^-1
end


function bilerp(x::Float64, y::Float64, Δx::Float64, Δy::Float64, x_samples::Vector{Float64}, y_samples::Vector{Float64}, LUT::Matrix{Float64})
    a, b = bracket(x, x_samples)
    p, q = bracket(y, y_samples)
    return (Δx * Δy)^-1 * (
    (x_samples[b] - x) * (y_samples[q] - y) * LUT[a, p] +
    (x_samples[b] - x) * (y - y_samples[p]) * LUT[a, q] + 
    (x - x_samples[a]) * (y_samples[q] - y) * LUT[b, p] +
    (x - x_samples[a]) * (y - y_samples[p]) * LUT[b, q])
end


# V = -μ / r_norm
# reference potential V_ref = -μ / R_E
# ΔV = V - V_ref
# geopotential height = ΔV / g_0
# some approximations are used:
# geodetic latitude is assumed to be launch site geodetic latitude
# truncated power series for latitude
function get_geopotential_height(r_norm::Float64)


    return μ * (-r_norm^-1 + R_E^-1) * g_0^-1
end


function get_soundspeed(T::Float64, molm::Float64)
    return sqrt(γ * 8.31446261815324 * T * molm^-1)
end


function get_normal_coefficient(Mach::Float64)
    local CN::Float64
    if Mach <= 0.0
        CN = 0.0
    else
        CN = lerp(Mach, ΔMach, Mach_samples, CN_LUT)
    end
    return CN
end


function get_aerodynamic_center(Mach::Float64)
    local CM::Vector{Float64}
    if Mach <= 0.0
        CM = [0.0; d_AVI - CM_LUT[1]; 0.0]
    else
        CM = [0.0; d_AVI - lerp(Mach, ΔMach, Mach_samples, CM_LUT); 0.0]
    end
    return CM
end


function get_axial_coefficient(Mach::Float64, α::Float64, engine_on::Bool)
    local CA::Float64
    if Mach <= 0.0
        CA = 0.0
    elseif engine_on == true
        CA = bilerp(Mach, α, ΔMach, Δα, Mach_samples, α_samples, CAon_LUT)
    else
        CA = bilerp(Mach, α, ΔMach, Δα, Mach_samples, α_samples, CAoff_LUT)
    end
    return CA
end


function get_inertia_matrix(m::Float64)
    local I::Matrix{Float64} = (I_dry * (m_wet_nom - m) + I_wet * (m - m_dry_nom)) * (m_wet_nom - m_dry_nom)^-1
    return I
end


function get_inertia_matrix_inv(m::Float64)
    local I_inv::Matrix{Float64} = (I_dry_inv * (m_wet_nom - m) + I_wet_inv * (m - m_dry_nom)) * (m_wet_nom - m_dry_nom)^-1
    return I_inv
end


function get_cg(m::Float64)
    local cg::Vector{Float64} = (cg_dry * (m_wet_nom - m) + cg_wet * (m - m_dry_nom)) * (m_wet_nom - m_dry_nom)^-1
    return cg
end


# primed vectors are body frame
function rail_function!(S::Vector{Float64}, t::Float64, output, params, i::Integer, j::Integer)
    local r_ENU::Vector{Float64} = S[1:3] # position vector
    local v_ENU::Vector{Float64} = S[4:6] # velocity vector
    # rotations are world -> body
    local q::Vector{Float64} = S[7:10] # quaternion
    local ω_ENU::Vector{Float64} = S[11:13] # body frame angular velocity
    local m_fuel::Float64 = S[14]
    local m_ox::Float64 = S[15]
    local q_NED::Vector{Float64} = S[16:19]


    local AM::Matrix{Float64} = q2AM(q)
    local DCM::Matrix{Float64} = q2DCM(q)
    
    
    local r_ECEF::Vector{Float64} = (pad_ENU2ECEF_mat * r_ENU) + pad_r_ECEF
    local v_ECEF::Vector{Float64} = (pad_ENU2ECEF_mat * v_ENU)
    local r_norm = norm(r_ECEF)
    local ω_Earth::Vector{Float64}
    local dωdt_Earth::Vector{Float64}
    ω_Earth, dωdt_Earth = get_ωα_ECEF(t[1], Δτ_iJ2k, Δτ_iY)
    local a_fictitious::Vector{Float64} = pad_ECEF2ENU_mat * (-(μ * r_norm^-2) * (r_ECEF * r_norm^-1)
    - cross(dωdt_Earth, r_ECEF)
    - 2.0 * cross(ω_Earth, v_ECEF)
    - cross(ω_Earth, cross(ω_Earth, r_ECEF)))
    local h::Float64 = r_ENU[3]
    local T::Float64
    local dTdh::Float64
    local molm::Float64
    local molmdh::Float64
    local ρ::Float64
    local dρdh::Float64
    local pstatic::Float64
    local dpdh::Float64
    T, dTdh = lerp_with_deriv(h, Δh, h_samples, temp_LUT)
    molm, molmdh = lerp_with_deriv(h, Δh, h_samples, molm_LUT)
    ρ, dρdh = lerp_with_deriv(h, Δh, h_samples, ρ_LUT)
    pstatic, dpdh = lerp_with_deriv(h, Δh, h_samples, pstatic_LUT)
    local soundspeed::Float64 = get_soundspeed(T, molm)
    local u′::Vector{Float64} = AM * (get_wind(φ_pad, λ_pad, r_ENU[1] + h_pad) - v_ENU) 
    local windspeed::Float64 = norm(u′)
    local Mach::Float64 = windspeed * soundspeed^-1
    local û′::Vector{Float64} = u′ * windspeed^-1
    local q_air::Float64 = 0.5 * ρ * windspeed^2
    local α::Float64 = acos(abs(û′[1]))
    local C_A::Float64 = get_axial_coefficient(Mach, α, true) # true means engine is on
    local C_N::Float64 = get_normal_coefficient(Mach)
    local û′_reduced::Vector{Float64} = û′ - [û′[1]; 0.0; 0.0]
    û′_reduced *= ((dot(û′_reduced, û′_reduced) > 0.0) ? norm(û′_reduced)^-1 : 0.0)
    local F_aero′::Vector{Float64} = q_air * S_ref * ((C_A * [sign(û′[1]); 0.0; 0.0]) + C_N * û′_reduced)
    local F_thrust′::Vector{Float64} = d̂_thrust * (thrust + S_nozzle * (pstatic - p_nozzle))
    local m::Float64 = m_dry_nom + m_fuel + m_ox + m_gas_nom


    local Σa_ENU::Vector{Float64} = a_fictitious + DCM * (m^-1 * (F_aero′ + F_thrust′))
    local cg::Vector{Float64} = get_cg(m)
    local C_M::Vector{Float64} = get_aerodynamic_center(Mach)
    local Inertia_inv::Matrix{Float64} = get_inertia_matrix_inv(m)
    local Σω̇_ENU ::Vector{Float64} = Inertia_inv * (cross(cg, (AM * a_fictitious)) + cross(C_M, F_aero′) + cross(r_thrust, F_thrust′))
    local q̇::Vector{Float64} = 0.5 * Ξ(q) * ω_ENU
    local ω_NED::Vector{Float64} = [ω_ENU[2]; ω_ENU[1]; -ω_ENU[3]]
    local q̇_NED::Vector{Float64} = 0.5 * Ξ(q_NED) * ω_NED
    local Ṡ::Vector{Float64} = [v_ENU; Σa_ENU; q̇; Σω̇_ENU; mdot_fuel; mdot_ox; q̇_NED]
    if j != 0
        local v_NED::Vector{Float64} = [v_ENU[2]; v_ENU[1]; -v_ENU[3]]
        local Σa_NED::Vector{Float64} = [Σa_ENU[2]; Σa_ENU[1]; -Σa_ENU[3]]
        local Σω̇_NED::Vector{Float64} = [Σω̇_ENU[2]; Σω̇_ENU[1]; -Σω̇_ENU[3]]
        local ω̇̂_readout::Vector{Float64} = [-Σω̇_NED[3]; Σω̇_NED[2]; Σω̇_NED[1]]
        local â_readout::Vector{Float64} = DCM * (m^-1 * (F_aero′ + F_thrust′))
        â_readout = [-â_readout[3]; â_readout[2]; â_readout[1]]
        output[1][1:15, j, i] = Ṡ[1:15]
        output[2][1:21, j, i] = [v_NED; Σa_NED; q̇_NED; Σω̇_NED; mdot_fuel; mdot_ox; ω̇̂_readout; â_readout]
        output[3][1:4, j, i] = v_ENU[3] * [dTdh; molmdh; dρdh; dpdh]
    end
    return Ṡ
end


# primed vectors are body frame
function engine_on_function!(S::Vector{Float64}, t::Float64, output, params, i::Integer, j::Integer)
    local r_ENU::Vector{Float64} = S[1:3] # position vector
    local v_ENU::Vector{Float64} = S[4:6] # velocity vector
    # rotations are world -> body
    local q::Vector{Float64} = S[7:10] # quaternion
    local ω_ENU::Vector{Float64} = S[11:13] # body frame angular velocity
    local m_fuel::Float64 = S[14]
    local m_ox::Float64 = S[15]
    local q_NED::Vector{Float64} = S[16:19]


    local AM::Matrix{Float64} = q2AM(q)
    local DCM::Matrix{Float64} = q2DCM(q)
    
    
    local r_ECEF::Vector{Float64} = (pad_ENU2ECEF_mat * r_ENU) + pad_r_ECEF
    local v_ECEF::Vector{Float64} = (pad_ENU2ECEF_mat * v_ENU)
    local r_norm = norm(r_ECEF)
    local ω_Earth::Vector{Float64}
    local dωdt_Earth::Vector{Float64}
    ω_Earth, dωdt_Earth = get_ωα_ECEF(t[1], Δτ_iJ2k, Δτ_iY)
    local a_fictitious::Vector{Float64} = pad_ECEF2ENU_mat * (-(μ * r_norm^-2) * (r_ECEF * r_norm^-1)
    - cross(dωdt_Earth, r_ECEF)
    - 2.0 * cross(ω_Earth, v_ECEF)
    - cross(ω_Earth, cross(ω_Earth, r_ECEF)))
    local h::Float64 = r_ENU[3]
    local T::Float64
    local dTdh::Float64
    local molm::Float64
    local molmdh::Float64
    local ρ::Float64
    local dρdh::Float64
    local pstatic::Float64
    local dpdh::Float64
    T, dTdh = lerp_with_deriv(h, Δh, h_samples, temp_LUT)
    molm, molmdh = lerp_with_deriv(h, Δh, h_samples, molm_LUT)
    ρ, dρdh = lerp_with_deriv(h, Δh, h_samples, ρ_LUT)
    pstatic, dpdh = lerp_with_deriv(h, Δh, h_samples, pstatic_LUT)
    local soundspeed::Float64 = get_soundspeed(T, molm)
    local u′::Vector{Float64} = AM * (get_wind(φ_pad, λ_pad, r_ENU[1] + h_pad) - v_ENU) 
    local windspeed::Float64 = norm(u′)
    local Mach::Float64 = windspeed * soundspeed^-1
    local û′::Vector{Float64} = u′ * windspeed^-1
    local q_air::Float64 = 0.5 * ρ * windspeed^2
    local α::Float64 = acos(abs(û′[1]))
    local C_A::Float64 = get_axial_coefficient(Mach, α, true) # true means engine is on
    local C_N::Float64 = get_normal_coefficient(Mach)
    local û′_reduced::Vector{Float64} = û′ - [û′[1]; 0.0; 0.0]
    û′_reduced *= ((dot(û′_reduced, û′_reduced) > 0.0) ? norm(û′_reduced)^-1 : 0.0)
    local F_aero′::Vector{Float64} = q_air * S_ref * ((C_A * [sign(û′[1]); 0.0; 0.0]) + C_N * û′_reduced)
    local F_thrust′::Vector{Float64} = d̂_thrust * (thrust + S_nozzle * (pstatic - p_nozzle))
    local m::Float64 = m_dry_nom + m_fuel + m_ox + m_gas_nom


    local Σa_ENU::Vector{Float64} = a_fictitious + DCM * (m^-1 * (F_aero′ + F_thrust′))
    local cg::Vector{Float64} = get_cg(m)
    local C_M::Vector{Float64} = get_aerodynamic_center(Mach)
    local Inertia_inv::Matrix{Float64} = get_inertia_matrix_inv(m)
    local Σω̇_ENU ::Vector{Float64} = Inertia_inv * (cross(cg, (AM * a_fictitious)) + cross(C_M, F_aero′) + cross(r_thrust, F_thrust′))
    local q̇::Vector{Float64} = 0.5 * Ξ(q) * ω_ENU
    local ω_NED::Vector{Float64} = [ω_ENU[2]; ω_ENU[1]; -ω_ENU[3]]
    local q̇_NED::Vector{Float64} = 0.5 * Ξ(q_NED) * ω_NED
    local Ṡ::Vector{Float64} = [v_ENU; Σa_ENU; q̇; Σω̇_ENU; mdot_fuel; mdot_ox; q̇_NED]
    if j != 0
        local v_NED::Vector{Float64} = [v_ENU[2]; v_ENU[1]; -v_ENU[3]]
        local Σa_NED::Vector{Float64} = [Σa_ENU[2]; Σa_ENU[1]; -Σa_ENU[3]]
        local Σω̇_NED::Vector{Float64} = [Σω̇_ENU[2]; Σω̇_ENU[1]; -Σω̇_ENU[3]]
        local ω̇̂_readout::Vector{Float64} = [-Σω̇_NED[3]; Σω̇_NED[2]; Σω̇_NED[1]]
        local â_readout::Vector{Float64} = DCM * (m^-1 * (F_aero′ + F_thrust′))
        â_readout = [-â_readout[3]; â_readout[2]; â_readout[1]]
        output[1][1:15, j, i] = Ṡ[1:15]
        output[2][1:21, j, i] = [v_NED; Σa_NED; q̇_NED; Σω̇_NED; mdot_fuel; mdot_ox; ω̇̂_readout; â_readout]
        output[3][1:4, j, i] = v_ENU[3] * [dTdh; molmdh; dρdh; dpdh]
    end
    return Ṡ
end







# ignore this
#=
# must be done in body frame
function get_body_aero_forces(u::Vector{Float64}, h::Float64, engine_on::Bool)
    local T::Float64
    local molm::Float64
    local ρ::Float64
    local p::Float64
    local soundspeed::Float64
    local airspeed::Float64
    local û::Vector{Float64}
    local û′::Vector{Float64}
    local Mach::Float64
    local α::Float64
    local CA::Float64
    local CN::Float64
    T = lerp(h, Δh, h_samples, temp_LUT)
    molm = lerp(h, Δh, h_samples, molm_LUT)
    ρ = lerp(h, Δh, h_samples, ρ_LUT)
    pstatic = lerp(h, Δh, h_samples, pstatic_LUT)
    soundspeed = get_soundspeed(T, molm)
    q = dot(u, u)
    airspeed = √q
    û = u * airspeed^-1
    Mach = airspeed * molm^-1
    q *= 0.5ρ
    α = acos(abs(û[1]))
    CA = get_axial_coefficient(Mach, α, engine_on)
    û′ = û - û[1]
    û′ *= ((dot(û′, û′) > 0.0) ? norm(û′)^-1 : 0.0)
    return q * S_ref * (CA * sign(û[1]) * [1.0; 0.0; 0.0] + CN * û′), get_aerodynamic_center(Mach)
end
=#

#=
function phase1_kinematics(r′::Vector{Float64}, v′::Vector{Float64},
    m::Float64, DCM::Matrix{Float64}, AM::Matrix{Float64}, engine_on::Bool,
    w′::Vector{Float64}, ω_Earth::Vector{Float64}, ω̇_Earth::Vector{Float64})
    r_norm = norm(r_ECEF)
    a_grav = -μ * r_norm^-2 * (r′ * r_norm^-1)
    F_aero, CM = DCM * get_body_aero_forces(AM * (w′ - v′), get_geopotential(r_norm), engine_on)
    F_thrust = [((lerp(h, Δh, h_samples, pstatic_LUT) - p_nozzle) * S_nozzle + thrust); 0.0; 0.0]
    Σa = a_grav + (m^-1 * (F_aero + F_thrust)) +
    (-cross(ω̇_Earth, r′)) + (-2.0 * cross(ω_Earth, v′)) + (-cross(ω_Earth, cross(ω_Earth, r′)))
    Στ = get_inertia_matrix_inv(m) * (m * cross(get_cg(m), a_grav) + cross(CM, F_aero) + cross(r_thrust, F_thrust))
    return Σa, Στ
end
=#


#=
function aero_forces(r, u_air, body_x̂, x_cg, flight_params)
    û_air = u_air / norm(u_air)
    Mach = norm(u_air) / 310.0
    q_atm = (1 // 2) * density(r) * sum(abs2.(u_air), dims = 1)[1]
    α = acos(abs(dot(body_x̂, -û_air)))
    aero_ẑ = cross(body_x̂, -û_air)
    aero_n̂ = cross(body_x̂, aero_ẑ)
    F_axial = axial_coefficient(α, Mach) * q_atm * A_ref * body_x̂
    F_normal = normal_coefficient(α, flight_params) * q_atm * A_ref * body_n̂
    F_aero = F_axial + F_normal
    C_m = pitching_coefficient(α, flight_params)
    x_ac = ((flight_params.x_cg / flight_params.x_ref) - (C_m / C_N)) * flight_params.x_ref
    τ_aero = cross((flight_params.x_cg - x_ac) * body_x̂, F_normal)
    return F_aero, τ_aero
end


function aero_kinematics(x̂, u)
    û = u * norm(u)^-1
    q = 0.5 * ρ(r) * dot(u, u)
    α = acos(abs(dot(x̂, -û)))


    q_atm = (1 // 2) * density(r) * sum(abs2.(u_air), dims = 1)[1]


end
=#



