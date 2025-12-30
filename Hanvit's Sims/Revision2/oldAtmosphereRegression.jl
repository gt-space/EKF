using LinearAlgebra
using NonlinearSolve


using SatelliteToolboxAtmosphericModels
using JLD2
include("Environment.jl")









function poly7_a(x, α)
    return α[1] + α[2] * x + α[3] * x^2 + α[4] * x^3 + α[5] * x^4 + α[6] * x^5 + α[7] * x^6 + α[8] * x^7
end


function poly7_b(x, β)
    β[1] * x + β[2] * x^2 + β[3] * x^3 + β[4] * x^4 + β[5] * x^5 + β[6] * x^6 + β[7] * x^7 + β[8] * x^8
end


function poly_objective_a(α, data)
    u, y = data
    ŷ = poly7_a.(u, (α, ))
    r = ŷ - y
    return r
end


function poly_objective_b(β, data)
    u, y = data
    ŷ = poly7_b.(u, (β, ))
    r = ŷ - y
    return r
end


function R_squared_a(u, y, C)
    RSS = sum(abs2.(y - exp.(poly7_a.(u, (C, )))), dims = 1)[1]
    TSS = sum(abs2.(y .- (sum(y, dims = 1)[1] * length(y)^-1)), dims = 1)[1]
    return 1.0 - RSS * TSS^-1
end


function R_squared_b(u, y, zero_point′, C)
    RSS = sum(abs2.(y - zero_point′ * exp.(poly7_b.(u, (C, )))), dims = 1)[1]
    TSS = sum(abs2.(y .- (sum(y, dims = 1)[1] * length(y)^-1)), dims = 1)[1]
    return 1.0 - RSS * TSS^-1
end


setprecision(128)


SpaceIndices.init()
high_samples::Vector{BigFloat} = collect(range(-50.0, 5.0 * 10.0^4, length = 2^8))
high_atm_data::Matrix{BigFloat} = collect(reduce(hcat, get_atmosphere_data.((35.34788, ), (-117.80683, ), Float64.(high_samples)))')


p_log::Vector{BigFloat} = log.(high_atm_data[:, 4])
zero_point::BigFloat = log(get_atmosphere_data(35.34788, -117.80683, 0.0)[4])
zero_point′::BigFloat = get_atmosphere_data(35.34788, -117.80683, 0.0)[4]
p_log_zeroed::Vector{BigFloat} = p_log .- zero_point


β_i::Vector{BigFloat} = BigFloat.(zeros(8))
problem_b = NonlinearLeastSquaresProblem(NonlinearFunction(poly_objective_b), β_i, (high_samples, p_log_zeroed))
solution_b = solve(problem_b)
β_f = solution_b.u


println(R_squared_b(high_samples, high_atm_data[:, 4], zero_point′, β_f))
println(sqrt(sum(abs2.((high_atm_data[:, 4] - zero_point′ * exp.(poly7_b.(high_samples, (β_f, )))) ./ high_atm_data[:, 4]), dims = 1)[1] .* length(p_log)^-1))


β_rounded = Float32.(β_f)
samples_rounded = Float32.(high_samples)
resultsf32::Vector{Float32} = exp.(Float32.(poly7_b.(samples_rounded, (β_rounded, )))) * Float32(zero_point′)


using GLMakie
fig = Figure()
ax = Axis(fig[1, 1], title = "Test Case", xlabel = "Altitude (metress)", ylabel = "Pressure (Pascal)")
plot_samples = collect(range(-10.0^3, 7.0 * 10.0^4, length = 2^8))
actual_values = collect(reduce(hcat, get_atmosphere_data.((35.34788, ), (-117.80683, ), plot_samples))')[:, 4]
lines!(ax, plot_samples, actual_values, label = "Raw Data")
lines!(ax, plot_samples, exp.(poly7_b.(plot_samples, (β_f, ))) * zero_point′, label = "Fit")
Legend(fig, ax)
fig


#=
α_i::Vector{BigFloat} = [log(10.0^5); zeros(7)]
problem_a = NonlinearLeastSquaresProblem(NonlinearFunction(poly_objective_a), α_i, (high_samples, p_log))
solution_a = solve(problem_a)
α_f = solution_a.u


println(R_squared_a(high_samples, high_atm_data[:, 4], α_f))
println(sqrt(sum(abs2.((exp.(p_log) - exp.(poly7_a.(high_samples, (α_f, )))) ./ exp.(p_log)), dims = 1)[1] .* length(p_log)^-1))


using GLMakie
fig = Figure()
ax = Axis(fig[1, 1], title = "Test Case", xlabel = "Altitude (meters)", ylabel = "Log Pressure")
lines!(ax, high_samples, p_log, label = "Raw Data")
lines!(ax, high_samples, (poly7_a.(high_samples, (α_f, ))), label = "Fit")
Legend(fig, ax)
fig
=#






#=
sample_means = 0.5 * (high_samples[2:end] + high_samples[1:end - 1])
differences = (log10.(high_atm_data[2:end, 4]) - log10.(high_atm_data[1:end - 1, 4])) ./ (high_samples[2:end] - high_samples[1:end - 1])


fig = Figure()
ax = Axis(fig[1, 1], title = "Pressure Derivative Plot", xlabel = "Altitude (metres)", ylabel = "Pressure Derivative (Pascal per metre)")
lines!(ax, sample_means, differences)
fig
=#



