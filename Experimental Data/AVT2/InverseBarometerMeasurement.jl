global α::Vector{Float32} = [
    5.0122185;
    -4.9929004e-5;
    -5.415637e-10;
    -3.837231e-14;
    2.55155e-18;
    -5.321706e-23;
    4.813401e-28;
    -1.6294356e-33
]
global A_p_ini::Float32 = (BigFloat(8.31446261815324) * BigFloat(291.355857442414)) * (BigFloat(-9.80665) * BigFloat(0.0289604202420187))^-1
global p_0_inv::Float32 = (BigFloat(1.22795595945348) * (BigFloat(8.31446261815324) / BigFloat(0.0289604202420187)) * BigFloat(291.355857442414))^-1
global A_p_ini64::Float64 = (BigFloat(8.31446261815324) * BigFloat(291.355857442414)) * (BigFloat(-9.80665) * BigFloat(0.0289604202420187))^-1
global p_0_inv64::Float64 = (BigFloat(1.22795595945348) * (BigFloat(8.31446261815324) / BigFloat(0.0289604202420187)) * BigFloat(291.355857442414))^-1


function poly7(x, α)
    return α[1] + α[2] * x + α[3] * x^2 + α[4] * x^3 + α[5] * x^4 + α[6] * x^5 + α[7] * x^6 + α[8] * x^7
end


function poly7deriv(x, α)
    return α[2] + 2.0 * α[3] * x + 3.0 * α[4] * x^2 + 4.0 * α[5] * x^3 + 5.0 * α[6] * x^4 + 6.0 * α[7] * x^5 + 7.0 * α[8] * x^6
end


function get_pressure(h::Float32, a::Vector{Float32})
    p::Float32 = 10.0^(a[1] + a[2] * h + a[3] * h^2 + a[4] * h^3 + a[5] * h^4 + a[6] * h^5 + a[7] * h^6 + a[8] * h^7)
    return p
end


function get_pressure′(h::Float32, a::Vector{Float32})
    p::Float32 = log(10.0) * (a[2] + 2.0 * a[3] * h + 3.0 * a[4] * h^2 + 4.0 * a[5] * h^3 + 5.0 * a[6] * h^4 + 6.0 * a[7] * h^5 + 7.0 * a[8] * h^6) * 10.0^(a[1] + a[2] * h + a[3] * h^2 + a[4] * h^3 + a[5] * h^4 + a[6] * h^5 + a[7] * h^6 + a[8] * h^7)
    return p
end


function invert_measurement_guess(p::Float32)
    return h::Float32 = A_p_ini * log(p * p_0_inv)
end
function invert_measurement_guess(p::Float64)
    return h::Float64 = A_p_ini64 * log(p * p_0_inv64)
end


# only works to 115 km
function p2alt(p::Float32)
    local h::Float32
    if (p < 105240.625f0) && (p > 100510.31f0)
        foo(h) = poly7(h, α) - log10(p)
        a::Float32 = -200.0
        b::Float32 = 200.0
        local c::Float32
        for i = 1:16
            c = 0.5 * (a + b)
            h = c
            if abs(b - a) < 10.0^-5
                break
            end
            if foo(c) == 0.0
                break
            elseif foo(c) > 0.0
                a = c
            elseif foo(c) < 0.0
                b = c
            end
        end
    else
    h = invert_measurement_guess(p)
    logp = log10(p)
        for i = 1:16
            Δh::Float32 = (poly7(h, α) - logp) * poly7deriv(h, α)^-1
            h -= Δh
            if abs(Δh) < 10.0^-5
                break
            end
        end
    end
    return h
end


# only works to 115 km
function p2alt(p::Float64)
    local h::Float64
    if (p < 105240.625) && (p > 100510.31)
        foo(h) = poly7(h, α) - log10(p)
        a::Float64 = -200.0
        b::Float64 = 200.0
        local c::Float64
        for i = 1:256
            c = 0.5 * (a + b)
            h = c
            if abs(b - a) < 10.0^-13
                break
            end
            if foo(c) == 0.0
                break
            elseif foo(c) > 0.0
                a = c
            elseif foo(c) < 0.0
                b = c
            end
        end
    else
    h = invert_measurement_guess(p)
    logp = log10(p)
        for i = 1:256
            Δh::Float64 = (poly7(h, α) - logp) * poly7deriv(h, α)^-1
            h -= Δh
            if abs(Δh) < 10.0^-13
                break
            end
        end
    end
    return h
end


#=
#test case
using BenchmarkTools
using GLMakie


h_sample = Float32.(collect(range(-2.0, 115.0, length = 512000)) * 10.0^3)
p_sample = get_pressure.(h_sample, (α, ))
predicted_h = p2alt.(p_sample)
Δh = h_sample - predicted_h
Δh_rel = Δh ./ h_sample


fig = Figure()
ax = Axis(fig[1, 1:2], title = "Barometer Function", xlabel = "height (metres)", ylabel = "Pressure (Pascal)")
lines!(ax, h_sample, p_sample, label = "True Height")
lines!(ax, predicted_h, p_sample, label = "Predicted Height")
Legend(fig, ax)


ax2 = Axis(fig[2, 1], title = "Absolute Altitude Estimation Error", xlabel = "Pressure", ylabel = "Absolute Error (metres)")
lines!(ax2, p_sample, Δh)
fig


ax3 = Axis(fig[2, 2], title = "Relative Altitude Estimation Error", xlabel = "Pressure", ylabel = "Relative Error")
lines!(ax3, p_sample, Δh_rel)
fig


@benchmark p2alt($(Float32(117.0rand() - 2.0)))
=#


