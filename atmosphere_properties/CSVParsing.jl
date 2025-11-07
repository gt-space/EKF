using CSV
using DataFrames
using Plots
using NonlinearSolve
using LinearAlgebra


function poly7(x, α)
    return α[1] + α[2] * x + α[3] * x^2 + α[4] * x^3 + α[5] * x^4 + α[6] * x^5 + α[7] * x^6 + α[8] * x^7
end


function exp10_poly7(h::Float32, α::Vector{Float32})
    y::Float32 = (Float32(10.0))^(poly7(h, α))
    return y
end


function poly_obj(α, data)
    h, y = data
    ŷ = poly7.(h, (α, )) # compute 7th order polynomial using the given 
    r = y - ŷ # residuals
    return r
end


function get_coeffs(α_i, h, y)
    problem = NonlinearLeastSquaresProblem(NonlinearFunction(poly_obj), α_i, (h, y))
    solution = solve(problem)
    α_f = solution.u
    return α_f
end


function μ_const_p(data)
    return sum(data, dims = 1)[1] * length(data)^-1
end


function R_squared(h, y, α)
    RSS = sum(abs2.(y - exp10_poly7.(Float32.(h), (Float32.(α), ))), dims = 1)[1]
    TSS = sum(abs2.(y .- μ_const_p(y)), dims = 1)[1]
    return 1.0 - RSS * TSS^-1
end


setprecision(128)


CSV_path = normpath(joinpath((@__FILE__), raw"..\pyatmos_output.csv"))
f = CSV.read(CSV_path, DataFrame)
h_mat = BigFloat.(collect(f[1, :])) * 10.0^3
ρ_mat = BigFloat.(collect(f[2, :]))
T_mat = BigFloat.(collect(f[3, :]))
M_mat = BigFloat.(collect(f[4, :]))
p_mat = ρ_mat .* (BigFloat(8.31446261815324) ./ M_mat) .* T_mat
ρ_0::BigFloat = 1.22795595945348
p_0 = ρ_0 * (BigFloat.(8.31446261815324) / BigFloat.(0.0289604202420187)) * BigFloat.(291.355857442414)
α_i = BigFloat.([
    log10(p_0);
    0.0;#1
    0.0;#2
    0.0;#3
    0.0;#4
    0.0;#5
    0.0;#6
    0.0#7
])
i_min = 1
i_max = 821


h̄_mat = h_mat[i_min:i_max]
ρ̄_mat = ρ_mat[i_min:i_max]
p̄_mat = p_mat[i_min:i_max]
ρ_log = log10.(ρ̄_mat)
p_log = log10.(p̄_mat)


α_p = get_coeffs(α_i, h̄_mat, p_log)
p_rsquared = R_squared(h̄_mat, p̄_mat, α_p)
println(p_rsquared)
α_p_f32 = Float32.(α_p)
println(α_p_f32)




