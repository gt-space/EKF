using CSV
using DataFrames
using Plots
using NonlinearSolve
using LinearAlgebra


CSV_path = normpath(joinpath((@__FILE__), raw"..\pyatmos_output.csv"))
f = CSV.read(CSV_path, DataFrame)
h_mat = collect(f[1, :]) * 10^3
ρ_mat = collect(f[2, :])
T_mat = collect(f[3, :])
M_mat = collect(f[4, :])
p_mat = ρ_mat .* (8.31446261815324 ./ M_mat) .* T_mat


function polynomial_regression(α, data)
    h, log_ρ = data
    prediction = polynomial_function.(h, (α, ))
    residuals = log_ρ - prediction
    return residuals
end


function polynomial_function(x, α)
    return (transpose(α) * [1.0; x; x^2; x^3; x^4; x^5; x^6; x^7; x^8; x^9; x^10; x^11; x^12; x^13; x^14; x^15])[1]
end


α_i = [
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0
]
ρ_0 = 1.22795595945348
p_0 = ρ_0 * (8.31446261815324 ./ 0.0289604202420187) .* 291.355857442414
max_index = 1621
small_h_mat = h_mat[1:max_index]
small_ρ_mat = ρ_mat[1:max_index]
small_p_mat = p_mat[1:max_index]


ρ_problem = NonlinearLeastSquaresProblem(NonlinearFunction(polynomial_regression), α_i, (small_h_mat, log10.(small_ρ_mat) .- log10.(ρ_0)))
ρ_solution = solve(ρ_problem)
ρ_coefficients = ρ_solution.u
p_problem = NonlinearLeastSquaresProblem(NonlinearFunction(polynomial_regression), α_i, (small_h_mat, log10.(small_p_mat) .- log10.(p_0)))
p_solution = solve(p_problem)
p_coefficients = p_solution.u


ρ_mean = sum(small_ρ_mat, dims = 1)[1] / max_index
ρ_RSS = sum(abs2.(small_ρ_mat .- ρ_0 * 10.0.^(polynomial_function.(small_h_mat, (ρ_coefficients, )))), dims = 1)[1]
ρ_TSS = sum(abs2.(small_ρ_mat .- ρ_mean), dims = 1)[1]
ρ_Rsq = 1.0 - ρ_RSS / ρ_TSS
p_mean = sum(small_p_mat, dims = 1)[1] / max_index
p_RSS = sum(abs2.(small_p_mat .- p_0 * 10.0.^(polynomial_function.(small_h_mat, (p_coefficients, )))), dims = 1)[1]
p_TSS = sum(abs2.(small_p_mat .- p_mean), dims = 1)[1]
p_Rsq = 1.0 - p_RSS / p_TSS


println(ρ_coefficients)
ρ_plot = plot(h_mat, log10.(ρ_mat) .- log10.(ρ_0), label = "NRLMSISE-00")
plot!(ρ_plot, h_mat, polynomial_function.(h_mat, (ρ_coefficients, )), label = "curve fit", xlims = (-15000, 200000), ylims = (-15, 0.2))


println(p_coefficients)
p_plot = plot(h_mat, log10.(p_mat) .- log10.(p_0), label = "NRLMSISE-00")
plot!(p_plot, h_mat, polynomial_function.(h_mat, (p_coefficients, )), label = "curve fit", xlims = (-10 * 10^3, 200 * 10^3), ylims = (-10, 0.5))



println(ρ_coefficients)
println(ρ_0)
println(ρ_Rsq)
println(p_coefficients)
println(p_0)
println(p_Rsq)



