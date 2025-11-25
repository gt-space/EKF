


using LinearAlgebra


include("Rotations.jl")
include("CoordinateTransformations.jl")
include("Kinematics.jl")
include("DESolver.jl")
include("NonlinearSolve.jl")
include("Environment.jl")


current_epoch = Epoch{UniversalTime}(2026, 1, 10, 12, 00, 0.0)
reference_epoch = UT1Epoch(0.0seconds, origin = :j2000)
year_epoch = Epoch{UniversalTime}(2000, 1, 1, 0, 0, 0.0)
global Δτ_iJ2k::Float64 = value(seconds(current_epoch - reference_epoch))
global Δτ_iY::Float64 = value(seconds(current_epoch - year_epoch))
get_ωα_ECEF(0.0, Δτ_iJ2k, Δτ_iY)


φ_pad::Float64 = 35.34788
λ_pad::Float64 = -117.80683
h_pad::Float64 = 625.0

global pad_ENU2ECEF_mat = [
    -sin(λ_pad) -sin(φ_pad)cos(λ_pad) cos(φ_pad)cos(λ_pad);
    cos(λ_pad)  -sin(φ_pad)sin(λ_pad) cos(φ_pad)sin(λ_pad);
    0.0          cos(φ_pad)            sin(φ_pad)
]
global pad_ECEF2ENU_mat = [
    -sin(λ_pad)           cos(λ_pad)            0.0;
    -sin(φ_pad)cos(λ_pad) -sin(φ_pad)sin(λ_pad) cos(φ_pad);
    cos(φ_pad)cos(λ_pad)  cos(φ_pad)sin(λ_pad)  sin(φ_pad)
]
global pad_r_ECEF = lla2ecef(φ_pad, λ_pad, h_pad)


l_rail::Float64 = 11.05


yaw_i::Float64   = 0.0
pitch_i::Float64 = 0.0
roll_i::Float64  = rand() * 2.0π
rail_hat::Vector{Float64} = zyx2DCM(yaw_i, pitch_i, roll_i)[:, 1]


r_0::Vector{Float64} = (l_OA - d_AVI) * rail_hat
S::Vector{Float64} = [r_0; zeros(3); zyx2q(yaw_i, pitch_i, roll_i); zeros(3); m_fuel_nom; m_ox_nom; zyx2q(yaw_i, pitch_i, roll_i)]
t::Vector{Float64} = [0.0]


n_guess::Int = 200
phase1_initial_states::Matrix{Float64} = zeros(19, n_guess)
phase1_initial_states[1:19, 1] = S
phase1_initial_times::Matrix{Float64} = zeros(1, n_guess)
phase1_initial_times[1] = 0.0
phase1_t_f::Float64 = 0.0
phase1_Δt::Float64 = 0.005
phase1_ε::Matrix{Float64} = zeros(19, n_guess)
phase1_output = (zeros(15, 19, n_guess), zeros(21, 19, n_guess), zeros(4, 19, n_guess))


global i = 1
while true
    phase1_ε[1:19, i] = step_sol_v9!(S, t, phase1_output, (), phase1_Δt, rail_function!, i)
    if norm(S[1:3] - r_0) > l_rail
        my_interpolant = phase1_output[1][1:3, 1:19, i]
        local foo(r::Vector{Float64}) = norm(r - r_0) - l_rail
        local foo′(u::Float64) = foo(interpolate_embed9(u, my_interpolant))
        u_terminate::Tuple{Float64, Float64} = ITP(0.0, 1.0, foo′; y_a = foo(phase1_initial_states[1:3, i]), y_b = foo(S[1:3]))
        phase1_t_f = phase1_initial_times[i] + phase1_Δt * ((1 // 2) * sum(u_terminate))
        break
    end
    i += 1
    phase1_initial_states[:, i] = S
    phase1_initial_times[i] = t[1]
end







# output is tuple
# each tuple entry is one type of data eg state deriv, fake state deriv
# each data type entry is n x 1 vector
# overall forms n x interpolant order x timestep count matrix
# 19 terms in vern9 interpolant 


# p, T, q, Mach, α
# force
# the other quaterion


# actual state
# fake state
# atmosphere properties












