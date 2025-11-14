



r_0 = lla2ecef(35.34788, -117.80683, 625.0)
v_0 = [0.0; 0.0; 0.0]
q_0 = zyx2q()
ω_0 = [0.0; 0.0; 0.0]
m_fuel_0 = m_fuel_nom
m_ox_0 = m_ox_nom


S = [r_0; v_0; q_0; ω_0; m_fuel_0; m_ox_0]


phase1_S_i::Matrix{Float64} = zeros()
phase1_t_f::Float64 = 0.0

while true
    step!(S, t, interpolant, params, Δt_1, f, i)

    if norm(S[1:3] - r_0) > l_rail
        my_interpolant = interpolant[1][i][1:3, :]
        break
    end
    i += 1
    phase1_S_i[:, i] = S
end













# output is tuple
# each tuple entry is one type of data eg state deriv, q, pressure
# each data type entry is n x 1 vector
# overall forms n x interpolant order x timestep count matrix
# 19 terms in vern9 interpolant 


# p, T, q, Mach, α
# force
# the other quaterion



function step!(S, t, output, params, Δt, f, i)
end












