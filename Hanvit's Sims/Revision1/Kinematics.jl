using LinearAlgebra


using JLD2


using GLMakie


perf_dict = load(normpath(joinpath((@__FILE__), raw"..\PerformanceData.jld2")))


global μ::Float64 = 3.986004418 * 10.0^14
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


global m_dry::Float64 = 128.48
global m_ox::Float64 = 80.85
global m_fuel::Float64 = 40.25
global m_gas::Float64 = 15.0
global m_wet::Float64 = m_dry + m_ox + m_fuel + m_gas
# simulation body frame y is positive towards nose
# z in CAD is x in simulation body frame
# y in CAD is z in simulation body frame
l_OA::Float64 = 7.23 + 0.15875
global d_AVI::Float64 = 41.03 * 0.0254
global r_thrust::Vector{Float64} = [0.0; -l_OA; 0.0] + [0.0; d_AVI; 0.0]
global cg_dry::Vector{Float64} = r_thrust + [0.00127; 4.18719; 0.000762]
global cg_wet::Vector{Float64} = r_thrust + [-0.000254; 3.681; 0.000508]
global I_dry::Matrix{Float64} = [
    3547.558 -2.4294 -0.01225;
    -2.4294 6.3856 -1.0125;
    -0.01225 -1.0125 3547.5412]
    + m_dry * (dot(cg_dry, cg_dry) * I(3) - cg_dry * transpose(cg_dry))
global I_wet::Matrix{Float64} = [
    4588.36 -2.194 -0.0119;
    -2.194 10.1713 -0.579;
    -0.0119 -0.579 4588.36]
    + m_wet * (dot(cg_wet, cg_wet) * I(3) - cg_wet * transpose(cg_wet))
global I_dry_inv = inv(I_dry)
global I_wet_inv = inv(I_wet) 
global S_ref::Float64 = 0.04047089
global thrust::Float64 = 11120.554
global S_nozzle::Float64 = 0.04047089
global p_nozzle::Float64 = 68947.0


function bracket(x::Float64, samples::Vector{Float64})
    local a::Int
    local b::Int
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


function bilerp(x::Float64, y::Float64, Δx::Float64, Δy::Float64, x_samples::Vector{Float64}, y_samples::Vector{Float64}, LUT::Matrix{Float64})
    a, b = bracket(x, x_samples)
    p, q = bracket(y, y_samples)
    return (Δx * Δy)^-1 * (
    (x_samples[b] - x) * (y_samples[q] - y) * LUT[a, p] +
    (x_samples[b] - x) * (y - y_samples[p]) * LUT[a, q] + 
    (x - x-samples[a]) * (y_samples[q] - y) * LUT[b, p] +
    (x - x_samples[a]) * (y - y_samples[p]) * LUT[b, q])
end


# V = -μ / r_norm
# reference potential V_ref = -μ / R_E
# ΔV = V - V_ref
# geopotential height = ΔV / g_0
function get_geopotential(r_norm::Float64)
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
    local CM::Float64
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
    return (I_dry * (m_wet - m) + I_wet * (m - m_dry)) * (m_wet - m_dry)^-1
end


function get_inertia_matrix_inv(m::Float64)
    return (I_dry_inv * (m_wet - m) + I_wet_inv * (m - m_dry)) * (m_wet - m_dry)^-1
end


function get_cg(m::Float64)
    return (cg_dry * (m_wet - m) + cg_wet * (m - m_dry)) * (m_wet - m_dry)^-1
end


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



