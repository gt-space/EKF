

"""
    bracket(u::Float64, sp::Vector{Float64})

Search the ascending list of sample points `sp` for the largest value smaller than `u` and the smallest value larger than `u`, and return their indices.

...
# Arguments
- `u::Float64`: The value to search the list with
- `sp::Vector{Float64}`: The sample points to be searched, in ascending order

...
"""
function bracket(u::Float64, sp::Vector{Float64})
    local a::Integer
    local b::Integer
    if u <= sp[1]
        a = 1
        b = 2
    elseif u >= sp[end]
        b = length(sp)
        a = b - 1
    else
        a = findlast(u .> sp)
        b = a + 1
    end
    return a, b
end


"""
    lerp(u::Float64, Δ_u::Float64, u_sp::Vector{Float64}, LUT::Vector{Float64})

Linearly interpolate the look-up table `LUT` at the value `u`.

...
# Arguments
- `u::Float64`: The input value to search the LUT with
- `Δ_u::Float64`: The spacing between LUT sample points
- `u_sp::Vector{Float64}`: The input value sample points, in ascending order
- `LUT::Vector{Float64}`: The values to be interpolated

...
"""
function lerp(u::Float64, Δ_u::Float64, u_sp::Vector{Float64}, LUT::Vector{Float64})
    local i::Integer
    local j::Integer
    i, j = bracket(u, u_sp)
    return (LUT[i] * (u_sp[j] - u) + LUT[j] * (u - u_sp[i])) * Δ_u^-1
end
function lerp(u::Float64, u_sp::Vector{Float64}, LUT::Vector{Float64})
    local i::Integer
    local j::Integer
    local Δ_u::Float64
    i, j = bracket(u, u_sp)
    Δ_u = u_sp[j] - u_sp[i]
    return (LUT[i] * (u_sp[j] - u) + LUT[j] * (u - u_sp[i])) * Δ_u^-1
end
function lerp(u::Float64, Δ_u::Float64, u_sp::Vector{Float64}, LUT::Matrix{Float64})
    local i::Integer
    local j::Integer
    i, j = bracket(u, u_sp)
    return (LUT[:, i] * (u_sp[j] - u) + LUT[:, j] * (u - u_sp[i])) * Δ_u^-1
end
function lerp(u::Float64, u_sp::Vector{Float64}, LUT::Matrix{Float64})
    local i::Integer
    local j::Integer
    local Δ_u::Float64
    i, j = bracket(u, u_sp)
    Δ_u = u_sp[j] - u_sp[i]
    return (LUT[:, i] * (u_sp[j] - u) + LUT[:, j] * (u - u_sp[i])) * Δ_u^-1
end


"""
    lerp_with_deriv(u::Float64, Δ_u::Float64, u_sp::Vector{Float64}, LUT::Vector{Float64})

Linearly interpolate the look-up table `LUT` at the value `u`, returning the derivative as well.

...
# Arguments
- `u::Float64`: The input value to search the LUT with
- `Δ_u::Float64`: The spacing between LUT sample points
- `u_sp::Vector{Float64}`: The input value sample points, in ascending order
- `LUT::Vector{Float64}`: The values to be interpolated

...
"""
function lerp_with_deriv(u::Float64, Δ_u::Float64, u_sp::Vector{Float64}, LUT::Vector{Float64})
    local i::Integer
    local j::Integer
    i, j = bracket(u, u_sp)
    return (LUT[i] * (u_sp[j] - u) + LUT[j] * (u - u_sp[i])) * Δ_u^-1, (u_sp[j] - u_sp[i]) * Δ_u^-1
end
function lerp_with_deriv(u::Float64, u_sp::Vector{Float64}, LUT::Vector{Float64})
    local i::Integer
    local j::Integer
    local Δ_u::Float64
    i, j = bracket(u, u_sp)
    Δ_u = u_sp[j] - u_sp[i]
    return (LUT[i] * (u_sp[j] - u) + LUT[j] * (u - u_sp[i])) * Δ_u^-1, (u_sp[j] - u_sp[i]) * Δ_u^-1
end


"""
    bilerp(u::Float64, v::Float64, Δu::Float64, Δv::Float64, u_sp::Vector{Float64}, v_sp::Vector{Float64}, LUT::Matrix{Float64})

Bilinearly interpolate the look-up table `LUT` at the value `x`.

...
# Arguments
- `u::Float64`: The u-dimension input value to search the LUT with
- `v::Float64`: The v-dimension input value to search the LUT with
- `Δ_u::Float64`: The spacing between LUT sample points in the u-dimension
- `Δ_v::Float64`: The spacing between LUT sample points in the v-dimension
- `u_sp::Vector{Float64}`: The u-dimension input value sample points, in ascending order
- `v_sp::Vector{Float64}`: The v-dimension input value sample points, in ascending order
- `LUT::Matrix{Float64}`: The values to be interpolated

...
"""
function bilerp(u::Float64, v::Float64, Δ_u::Float64, Δ_v::Float64, u_sp::Vector{Float64}, v_sp::Vector{Float64}, LUT::Matrix{Float64})
    local i::Integer
    local j::Integer
    local k::Integer
    local l::Integer
    i, j = bracket(u, u_sp)
    k, l = bracket(v, v_sp)
    return (Δ_u * Δ_v)^-1 * (
    (u_sp[j] - u) * (v_sp[l] - v) * LUT[i, k] +
    (u_sp[j] - u) * (v - v_sp[k]) * LUT[i, l] + 
    (u - u_sp[i]) * (v_sp[l] - v) * LUT[j, k] +
    (u - u_sp[i]) * (v - v_sp[k]) * LUT[j, l])
end
function bilerp(u::Float64, v::Float64, u_sp::Vector{Float64}, v_sp::Vector{Float64}, LUT::Matrix{Float64})
    local i::Integer
    local j::Integer
    local k::Integer
    local l::Integer
    local Δ_u::Float64
    local Δ_v::Float64
    i, j = bracket(u, u_sp)
    k, l = bracket(v, v_sp)
    Δ_u = u_sp[j] - u_sp[i]
    Δ_v = v_sp[l] - v_sp[k]
    return (Δ_u * Δ_v)^-1 * (
    (u_sp[j] - u) * (v_sp[l] - v) * LUT[i, k] +
    (u_sp[j] - u) * (v - v_sp[k]) * LUT[i, l] + 
    (u - u_sp[i]) * (v_sp[l] - v) * LUT[j, k] +
    (u - u_sp[i]) * (v - v_sp[k]) * LUT[j, l])
end





