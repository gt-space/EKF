using LinearAlgebra


# scalar last formalism
# quaternion is double precision length 4 vector
# AM is world -> body
# DCM is body -> world
# DCM = AM^-1 = AM^T
# ω_I = AM^T ω = DCM ω
# ω is body frame angular velocity
# ω_I is world angular velocity


function get_cross_eqv(x::Vector{Float64})
    local cross_eqv_mat::Matrix{Float64} = [
        0.0 -x[3] x[2];
        x[3] 0.0 -x[1];
        -x[2] x[1] 0.0
    ]
    return cross_eqv_mat
end


function conj_q(q::Vector{Float64})
    local q′::Vector{Float64} = [-q[1]; -q[2]; -q[3]; q[4]]
    return q′
end


function inv_q(q::Vector{Float64})
    local q_inv::Vector{Float64} = conj_q(q) * dot(q, q)^-1
    return q_inv
end


function ⊗(q::Vector{Float64})
    local other_prod_mat::Matrix{Float64} = [
        q[4] q[3] -q[2] q[1];
        -q[3] q[4] q[1] q[2];
        q[2] -q[1] q[4] q[3];
        -q[1] -q[2] -q[3] q[4]
    ]
    return other_prod_mat
end
function ⊗(q̄::Vector{Float64}, q::Vector{Float64})
    local q′::Vector{Float64} = [
        q[4] * q̄[1:3] + q̄ * q[1:3] - cross(q̄[1:3], q[1:3]);
        q̄[4] * q[4] - dot(q̄[1:3], q[1:3])
    ]
    return q′
end


function ⊙(q::Vector{Float64})
    local Hamilton_prod_mat::Matrix{Float64} = [
        q[4] -q[3] q[2] q[1];
        q[3] q[4] -q[1] q[2];
        -q[2] q[1] q[4] q[3];
        -q[1] -q[2] -q[3] q[4]
    ]
    return Hamilton_prod_mat
end
function ⊙(q̄::Vector{Float64}, q::Vector{Float64})
    local q′::Vector{Float64} = [
        q[4] * q̄[1:3] + q̄ * q[1:3] + cross(q̄[1:3], q[1:3]);
        q̄[4] * q[4] - dot(q̄[1:3], q[1:3])
    ]
    return q′
end


function Ψ(q::Vector{Float64})
    local Ψ::Matrix{Float64} = [
        q[4] q[3] -q[2];
        -q[3] q[4] q[1];
        q[2] -q[1] q[4];
        -q[1] -q[2] -q[3]
    ]
    return Ψ
end


function Ξ(q::Vector{Float64})
    local Ξ::Matrix{Float64} = [
        q[4] -q[3] q[2];
        q[3] q[4] -q[1];
        -q[2] q[1] q[4];
        -q[1] -q[2] -q[3]
    ]
    return Ξ
end


function q2AM(q::Vector{Float64})
    local AM::Matrix{Float64} = Ξ(q)' * Ψ(q)
    return AM
end


function q2DCM(q::Vector{Float64})
    local DCM::Matrix{Float64} = q2AM(inv_q(q))
    return DCM
end


function zyx2q(ψ::Float64, θ::Float64, φ::Float64)
    local q::Vector{Float64} = [
        cos(0.5 * ψ) * cos(0.5 * θ) * sin(0.5 * φ) - sin(0.5 * ψ) * sin(0.5 * θ) * cos(0.5 * φ);
        cos(0.5 * ψ) * sin(0.5 * θ) * cos(0.5 * φ) + sin(0.5 * ψ) * cos(0.5 * θ) * sin(0.5 * φ);
        sin(0.5 * ψ) * cos(0.5 * θ) * cos(0.5 * φ) - cos(0.5 * ψ) * sin(0.5 * θ) * sin(0.5 * φ);
        cos(0.5 * ψ) * cos(0.5 * θ) * cos(0.5 * φ) + sin(0.5 * ψ) * sin(0.5 * θ) * sin(0.5 * φ)
    ]
    return q
end
function zyx2q(angles::Tuple{Float64, Float64, Float64})
    local ψ::Float64 = angles[1]
    local θ::Float64 = angles[2]
    local φ::Float64 = angles[3]
    local q::Vector{Float64} = [
        cos(0.5 * ψ) * cos(0.5 * θ) * sin(0.5 * φ) - sin(0.5 * ψ) * sin(0.5 * θ) * cos(0.5 * φ);
        cos(0.5 * ψ) * sin(0.5 * θ) * cos(0.5 * φ) + sin(0.5 * ψ) * cos(0.5 * θ) * sin(0.5 * φ);
        sin(0.5 * ψ) * cos(0.5 * θ) * cos(0.5 * φ) - cos(0.5 * ψ) * sin(0.5 * θ) * sin(0.5 * φ);
        cos(0.5 * ψ) * cos(0.5 * θ) * cos(0.5 * φ) + sin(0.5 * ψ) * sin(0.5 * θ) * sin(0.5 * φ)
    ]
    return q
end


function zyx2DCM(ψ::Float64, θ::Float64, φ::Float64)
    DCM::Matrix{Float64} = [
        (cos(ψ)cos(θ)) (cos(ψ)sin(θ)sin(φ) - sin(ψ)cos(φ)) (cos(ψ)sin(θ)cos(φ) + sin(ψ)sin(φ));
        (sin(ψ)cos(θ)) (sin(ψ)sin(θ)sin(φ) + cos(ψ)cos(φ)) (sin(ψ)sin(θ)cos(φ) - cos(ψ)sin(φ));
        -sin(θ) cos(θ)sin(φ) cos(θ)cos(φ)
    ]
    return DCM
end
function zyx2DCM(angles::Tuple{Float64, Float64, Float64})
    local ψ::Float64 = angles[1]
    local θ::Float64 = angles[2]
    local φ::Float64 = angles[3]
    DCM::Matrix{Float64} = [
        (cos(ψ)cos(θ)) (cos(ψ)sin(θ)sin(φ) - sin(ψ)cos(φ)) (cos(ψ)sin(θ)cos(φ) + sin(ψ)sin(φ));
        (sin(ψ)cos(θ)) (sin(ψ)sin(θ)sin(φ) + cos(ψ)cos(φ)) (sin(ψ)sin(θ)cos(φ) - cos(ψ)sin(φ));
        -sin(θ) cos(θ)sin(φ) cos(θ)cos(φ)
    ]
    return DCM
end


#=
function DCM2zyx(DCM::Matrix{Float64})
    local ψ::Float64 = atan(DCM[2, 1], DCM[1, 1])                       # yaw angle
    local θ::Float64 = atan(-DCM[3, 1], sqrt(1.0 - abs2(DCM[3, 1])))    # pitch angle
    local φ::Float64 = atan(DCM[3, 2], DCM[3, 3])                       # roll angle
    return ψ, θ, φ
end


function AM2q(AM::Matrix{Float64})
    local q::Vector{Float64} = q_inv(zyx2q(DCM2zyx(AM)))
    return q
end


function DCM2q(DCM::Matrix{Float64})
    local q::Vector{Float64} = zyx2q(DCM2zyx(DCM))
    return q
end


function q2zyx(q::Vector{Float64})
    local ψ::Float64
    local θ::Float64
    local φ::Float64
    ψ, θ, φ = DCM2zyx(q2DCM(q))
    return ψ, θ, φ
end
=#




