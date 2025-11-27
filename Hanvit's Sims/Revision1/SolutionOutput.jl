function find_start(u::Float64, samples::Vector{Float64})
    local a::Integer
    if u <= samples[1]
        a = 1
    elseif u >= samples[end]
        a = length(samples)
    else
        a = findlast(u .> samples)
    end
    return a
end


function interpolate_sol(u::Float64, Δu::Float64, sol_samples::Vector{Float64}, sol_states::Matrix{Float64}, sol_derivs::Array{Float64, 3}; u_f::Float64 = Inf64)
    if u > u_f
        error("out of bounds")
    end
    i_ini::Integer = find_start(u, sol_samples)
    u′::Float64 = (u - sol_samples[i_ini]) * Δu^-1
    interpolant::Vector{Float64} = sol_states[1:15, i_ini] + Δu * interpolate_embed9(u′, sol_derivs[:, 1:19, i_ini])
    return interpolant
end


function interpolate_sol2(u::Float64, Δu::Float64, sol_samples::Vector{Float64}, sol_states::Matrix{Float64}, sol_derivs::Array{Float64, 3}; u_f::Float64 = Inf64)
    if u > u_f
        error("out of bounds")
    end
    i_ini::Integer = find_start(u, sol_samples)
    u′::Float64 = (u - sol_samples[i_ini]) * Δu^-1
    return u′
end






