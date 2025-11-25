using BenchmarkTools

#= 
function benchmark function(x::Float64, p1::Real, p2::Real)
    A = log(abs(x - p2))
    B = (1 // 3) + sign(x) * abs(x)^(1 // 3) + x^3
    y = x - p1
    C = (0.5y)^2 + ceil(0.5y) - 0.5
    return -(A + B + C)
end
=#



# must call with y_a and y_b having different signs
function ITP(a::Float64, b::Float64, f::Function; y_a::Float64 = f(a), y_b::Float64 = f(b), n::Integer = 128, ε::Float64 = 10.0^-13, κ_1::Float64 = 0.1, κ_2::Float64 = 0.98 * (1.0 + Base.MathConstants.golden))
    invert::Bool = false
    # assumes the function f(x) is increasing, if it's decreasing we can just find the zeroes of -f(x)
    if y_a > y_b
        invert = true
        y_a *= -1.0
        y_b *= -1.0
    end
    n += ceil(log2(0.5(b - a) * ε^-1))
    i::Integer = 0
    local y::Float64
    local t::Float64
    local t_bi::Float64
    local t_rf::Float64
    local Δ_t::Float64
    local t_t::Float64
    for i = 0:n
        # bisection point calculation:
        t_bi = 0.5(a + b)
        # regula-falsi point calculation:
        t_rf = (y_b * a - y_a * b) * (y_b - y_a)^-1
        # truncation step:
        Δ_t = t_bi - t_rf
        σ = sign(Δ_t)
        δ = κ_1 * (b - a)^κ_2
        if δ > abs(Δ_t)
            t_t = t_bi
        else
            t_t = t_rf + σ * δ
        end
        # projection step:
        r_k = ε * 2.0^(n - i) - 0.5(b - a)
        if abs(t_t - t_bi) > r_k
            t = t_bi - σ * r_k
        else
            t = t_t
        end
        # update step:
        y = f(t)
        if invert == true
            y *= -1.0
        end
        if y > 0
            b = t
            y_b = y
        elseif y < 0
            a = t
            y_a = y
        else
            a = t
            b = t
        end
        i += 1
        if abs(b - a) < 2.0ε
            break
        end
    end
    return a, b
end
function ITP(a::Float64, b::Float64, f::Function, p; y_a::Float64 = f(a, p), y_b::Float64 = f(b, p), n::Integer = 128, ε::Float64 = 10.0^-13, κ_1::Float64 = 0.1, κ_2::Float64 = 0.98 * (1.0 + Base.MathConstants.golden))
    invert::Bool = false
    if y_a > y_b
        invert = true
        y_a *= -1.0
        y_b *= -1.0
    end
    n += ceil(log2(0.5(b - a) * ε^-1))
    i::Integer = 0
    local y::Float64
    local t::Float64
    local t_bi::Float64
    local t_rf::Float64
    local Δ_t::Float64
    local t_t::Float64
    for i = 0:n
        # bisection point calculation:
        t_bi = 0.5(a + b)
        # regula-falsi point calculation:
        t_rf = (y_b * a - y_a * b) * (y_b - y_a)^-1
        # truncation step:
        Δ_t = t_bi - t_rf
        σ = sign(Δ_t)
        δ = κ_1 * (b - a)^κ_2
        if δ > abs(Δ_t)
            t_t = t_bi
        else
            t_t = t_rf + σ * δ
        end
        # projection step:
        r_k = ε * 2.0^(n - i) - 0.5(b - a)
        if abs(t_t - t_bi) > r_k
            t = t_bi - σ * r_k
        else
            t = t_t
        end
        # update step:
        y = f(t, p)
        if invert == true
            y *= -1.0
        end
        if y > 0
            b = t
            y_b = y
        elseif y < 0
            a = t
            y_a = y
        else
            a = t
            b = t
        end
        i += 1
        if abs(b - a) < 2.0ε
            break
        end
    end
    return a, b
end


#=
function Auto_Newton(t::T, f::Function, p; l_s, n, ε) where {T<:Real}

end
=#



