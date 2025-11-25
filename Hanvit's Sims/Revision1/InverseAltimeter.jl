using BenchmarkTools


function get_pressure(h::Float32, a::Vector{Float32})
    p::Float32 = 10.0^(a[1] + a[2] * h + a[3] * h^2 + a[4] * h^3 + a[5] * h^4 + a[6] * h^5 + a[7] * h^6 + a[8] * h^7)
    return p
end


function get_pressure′(h::Float32, a::Vector{Float32})
    p::Float32 = log(10.0) * (a[2] + 2.0 * a[3] * h + 3.0 * a[4] * h^2 + 4.0 * a[5] * h^3 + 5.0 * a[6] * h^4 + 6.0 * a[7] * h^5 + 7.0 * a[8] * h^6) * 10.0^(a[1] + a[2] * h + a[3] * h^2 + a[4] * h^3 + a[5] * h^4 + a[6] * h^5 + a[7] * h^6 + a[8] * h^7)
    return p
end


function invert_measurement_guess(p::Float32)
    R::Float32 = 8.31446261815324
    ρ::Float32 = 1.22795595945348
    molm::Float32 = 0.0289604202420187
    T::Float32 = 291.355857442414
    A::Float32 = (R * T) * (-9.80665 * molm)^-1
    p_0_inv::Float32 = (BigFloat(1.22795595945348) * (BigFloat(8.31446261815324) / BigFloat(0.0289604202420187)) * BigFloat(291.355857442414))^-1
    return h::Float32 = A * log(p * p_0_inv) 
end


function p2alt(p::Float32)
    α::Vector{Float32} = [
        5.0122185;
        -4.9929004e-5;
        -5.415637e-10;
        -3.837231e-14;
        2.55155e-18;
        -5.321706e-23;
        4.813401e-28;
        -1.6294356e-33
    ]
    h::Float32 = invert_measurement_guess(p)
    for i = 1:100
        Δh::Float32 = (get_pressure(h, α) - p) * get_pressure′(h, α)^-1
        h -= Δh
        if abs(Δh) <= 10.0^-5
            break
        end
    end
    return h
end

#=
p::Float32 = Float32(rand()*100000.0)
h::Float32 = p2alt(p)
a::Vector{Float32} = [
        5.0122185;
        -4.9929004e-5;
        -5.415637e-10;
        -3.837231e-14;
        2.55155e-18;
        -5.321706e-23;
        4.813401e-28;
        -1.6294356e-33
    ]
println(p)
println(h)
println(get_pressure(h, a))
=#
#@benchmark p2alt($(Float32(rand()*100000.0)))

#get_pressure(invert_measurement_guess(p_0), a)







