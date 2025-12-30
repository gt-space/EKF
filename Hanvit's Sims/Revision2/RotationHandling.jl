function LinearAlgebra.cross(x::Vector{Float64})
    local cross_eqv_mat::Matrix{Float64} = [
        0.0 -x[3] x[2];
        x[3] 0.0 -x[1];
        -x[2] x[1] 0.0
    ]
    return cross_eqv_mat
end


"""
    conj_q(q::Vector{Float64})

Computes the conjugate of the input quaternion given as a vector in the scalar-last foramlism.

...
# Arguments
- `q::Vector{Float64}`: The scalar-last quaternion to be conjugated

...
"""
function conj_q(q::Vector{Float64})
    return q′::Vector{Float64} = [-q[1]; -q[2]; -q[3]; q[4]]
end


"""
    inv_q(q::Vector{Float64})

Computes the inverse of the input quaternion given as a vector in the scalar-last foramlism.

...
# Arguments
- `q::Vector{Float64}`: The scalar-last quaternion to be inverted

...
"""
function inv_q(q::Vector{Float64})
    return q_inv::Vector{Float64} = conj_q(q) * dot(q, q)^-1
end


"""
    ⊗(q̄::Vector{Float64}, q::Vector{Float64})

Return the product `q̄ ⊗ q` of the quaterions `q̄` and `q`. Call on a single quaternion `q` to return the equivalent matrix, `[q⊗] = [Ψ(q) q]`.

Not to be confused with the less useful but more notorious Hamilton product, `q̄ ⊙ q`.

...
"""
function ⊗(q̄::Vector{Float64}, q::Vector{Float64})
    return q′::Vector{Float64} = [
        q[4] * q̄[1:3] + q̄ * q[1:3] - cross(q̄[1:3], q[1:3]);
        q̄[4] * q[4] - dot(q̄[1:3], q[1:3])
    ]
end
⊗(q::Vector{Float64}) = [
    q[4] q[3] -q[2] q[1];
    -q[3] q[4] q[1] q[2];
    q[2] -q[1] q[4] q[3];
    -q[1] -q[2] -q[3] q[4]
]


"""
    ⊙(q̄::Vector{Float64}, q::Vector{Float64})

Return the Hamilton product `q̄ ⊙ q` of the quaterions `q̄` and `q`. Call on a single quaternion `q` to return the Hamilton product equivalent matrix, `[q⊙] = [Ξ(q) q]`.

Not to be confused with the more useful but less notorious product, `q̄ ⊗ q`.

...
"""
function ⊙(q̄::Vector{Float64}, q::Vector{Float64})
    return q′::Vector{Float64} = [
        q[4] * q̄[1:3] + q̄ * q[1:3] + cross(q̄[1:3], q[1:3]);
        q̄[4] * q[4] - dot(q̄[1:3], q[1:3])
    ]
end
⊙(q::Vector{Float64}) = [
    q[4] -q[3] q[2] q[1];
    q[3] q[4] -q[1] q[2];
    -q[2] q[1] q[4] q[3];
    -q[1] -q[2] -q[3] q[4]
]


"""
    Ψ(q::Vector{Float64}) 

Return the first three columns of `[q⊗]`.

...
"""
function Ψ(q::Vector{Float64})
    return Ψ::Matrix{Float64} = [
        q[4] q[3] -q[2];
        -q[3] q[4] q[1];
        q[2] -q[1] q[4];
        -q[1] -q[2] -q[3]
    ]
end


"""
    Ξ(q::Vector{Float64}) 

Return the first three columns of `[q⊙]`.

...
"""
function Ξ(q::Vector{Float64})
    return Ξ::Matrix{Float64} = [
        q[4] -q[3] q[2];
        q[3] q[4] -q[1];
        -q[2] q[1] q[4];
        -q[1] -q[2] -q[3]
    ]
end


"""
    q2AM(q::Vector{Float64})

Return an attitude matrix (AM) in the right-handed z, y, x representation matching the input quaternion.

...
"""
function q2AM(q::Vector{Float64})
    return AM::Matrix{Float64} = Ξ(q)' * Ψ(q)
end


"""
    q2DCM(q::Vector{Float64})

Return a direction cosine matrix (DCM) in the right-handed z, y, x representation matching the input quaternion.

...
"""
function q2DCM(q::Vector{Float64})
    return DCM::Matrix{Float64} = Ψ(q)' * Ξ(q)
end


"""
    zyx2q(ψ::Float64, θ::Float64, φ::Float64)

Return a quaternion `q` matching the Euler angles given in right-handed z, y, x Euler angle representation.

...
# Arguments
- `ψ::Float64`: Yaw angle (z-rotation)
- `θ::Float64`: Pitch angle (y′-rotation)
- `φ::Float64`: Roll angle (x′-rotation)

...
"""
function zyx2q(ψ::Float64, θ::Float64, φ::Float64)
    return q::Vector{Float64} = [
        cos(0.5 * ψ) * cos(0.5 * θ) * sin(0.5 * φ) - sin(0.5 * ψ) * sin(0.5 * θ) * cos(0.5 * φ);
        cos(0.5 * ψ) * sin(0.5 * θ) * cos(0.5 * φ) + sin(0.5 * ψ) * cos(0.5 * θ) * sin(0.5 * φ);
        sin(0.5 * ψ) * cos(0.5 * θ) * cos(0.5 * φ) - cos(0.5 * ψ) * sin(0.5 * θ) * sin(0.5 * φ);
        cos(0.5 * ψ) * cos(0.5 * θ) * cos(0.5 * φ) + sin(0.5 * ψ) * sin(0.5 * θ) * sin(0.5 * φ)
    ]
end
function zyx2q(angles::Tuple{Float64, Float64, Float64})
    local ψ::Float64 = angles[1]
    local θ::Float64 = angles[2]
    local φ::Float64 = angles[3]
    return q::Vector{Float64} = [
        cos(0.5 * ψ) * cos(0.5 * θ) * sin(0.5 * φ) - sin(0.5 * ψ) * sin(0.5 * θ) * cos(0.5 * φ);
        cos(0.5 * ψ) * sin(0.5 * θ) * cos(0.5 * φ) + sin(0.5 * ψ) * cos(0.5 * θ) * sin(0.5 * φ);
        sin(0.5 * ψ) * cos(0.5 * θ) * cos(0.5 * φ) - cos(0.5 * ψ) * sin(0.5 * θ) * sin(0.5 * φ);
        cos(0.5 * ψ) * cos(0.5 * θ) * cos(0.5 * φ) + sin(0.5 * ψ) * sin(0.5 * θ) * sin(0.5 * φ)
    ]
end


"""
    zyx2DCM(ψ::Float64, θ::Float64, φ::Float64)

Return a direction cosine matrix (DCM) matching the Euler angles given in right-handed z, y, x Euler angle representation.

...
# Arguments
- `ψ::Float64`: Yaw angle (z-rotation)
- `θ::Float64`: Pitch angle (y′-rotation)
- `φ::Float64`: Roll angle (x′-rotation)

...
"""
function zyx2DCM(ψ::Float64, θ::Float64, φ::Float64)
    return DCM::Matrix{Float64} = [
        (cos(ψ)cos(θ)) (cos(ψ)sin(θ)sin(φ) - sin(ψ)cos(φ)) (cos(ψ)sin(θ)cos(φ) + sin(ψ)sin(φ));
        (sin(ψ)cos(θ)) (sin(ψ)sin(θ)sin(φ) + cos(ψ)cos(φ)) (sin(ψ)sin(θ)cos(φ) - cos(ψ)sin(φ));
        -sin(θ) cos(θ)sin(φ) cos(θ)cos(φ)
    ]
end
function zyx2DCM(angles::Tuple{Float64, Float64, Float64})
    local ψ::Float64 = angles[1]
    local θ::Float64 = angles[2]
    local φ::Float64 = angles[3]
    return DCM::Matrix{Float64} = [
        (cos(ψ)cos(θ)) (cos(ψ)sin(θ)sin(φ) - sin(ψ)cos(φ)) (cos(ψ)sin(θ)cos(φ) + sin(ψ)sin(φ));
        (sin(ψ)cos(θ)) (sin(ψ)sin(θ)sin(φ) + cos(ψ)cos(φ)) (sin(ψ)sin(θ)cos(φ) - cos(ψ)sin(φ));
        -sin(θ) cos(θ)sin(φ) cos(θ)cos(φ)
    ]
end




