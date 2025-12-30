using LinearAlgebra
using NonlinearSolve


using JLD2
using GLMakie
using Statistics


include("InverseBarometerMeasurement.jl")


function constant_fit_function(α::Vector{Float64}, data)
    x, y = data
    r = y .- α[1]
    return r
end


function linear_fit_function(α::Vector{Float64}, data)
    x, y = data
    r = (α[1] * x .+ α[2]) - y
    return r
end


function get_constant_fit(x::Vector{Float64}, y::Vector{Float64})
    my_problem = NonlinearLeastSquaresProblem(NonlinearFunction(constant_fit_function), [0.0], (x, y))
    my_solution = solve(my_problem)
    μ = (my_solution.u)[1]
    return μ
end


function get_linear_fit(x::Vector{Float64}, y::Vector{Float64})
    my_problem = NonlinearLeastSquaresProblem(NonlinearFunction(linear_fit_function), [0.0; 0.0], (x, y))
    my_solution = solve(my_problem)
    lin_fit = my_solution.u
    return lin_fit
end


function get_mean(x::Vector{Float64}, P::Vector{Float64})
    return dot(x, P)
end
function get_mean(x::Vector{Float64})
    return sum(x) * length(x)^-1
end


function get_covariance(x::Vector{Float64}, y::Vector{Float64}, P::Vector{Float64}, μ_x::Float64, μ_y::Float64)
    return get_mean(P, (x .- μ_x) .* (y .- μ_y))
end
function get_covariance(x::Vector{Float64}, y::Vector{Float64}, P::Vector{Float64})
    local μ_x::Float64 = get_mean(x, P)
    local μ_y::Float64 = get_mean(y, P)
    return get_mean(P, (x .- μ_x) .* (y .- μ_y))
end
function get_covariance(x::Vector{Float64}, y::Vector{Float64}, μ_x::Float64, μ_y::Float64)
    return dot((x .- μ_x), (y .- μ_y)) * length(x)^-1
end
function get_covariance(x::Vector{Float64}, y::Vector{Float64})
    local μ_x::Float64 = get_mean(x)
    local μ_y::Float64 = get_mean(y)
    return dot((x .- μ_x), (y .- μ_y)) * length(x)^-1
end


function get_covariance_mat(v::Matrix{Float64}, P::Vector{Float64}, μ_v::Vector{Float64})
    xx = get_covariance(v[1, :], v[1, :], P, μ_v[1], μ_v[1])
    yy = get_covariance(v[2, :], v[2, :], P, μ_v[2], μ_v[2])
    zz = get_covariance(v[3, :], v[3, :], P, μ_v[3], μ_v[3])
    xy = get_covariance(v[2, :], v[1, :], P, μ_v[2], μ_v[1])
    xz = get_covariance(v[3, :], v[1, :], P, μ_v[3], μ_v[1])
    yz = get_covariance(v[3, :], v[2, :], P, μ_v[3], μ_v[2])
    return [
        xx xy xz;
        xy yy yz;
        xz yz zz
    ]
end
function get_covariance_mat(v::Matrix{Float64}, μ_v::Vector{Float64})
    xx = get_covariance(v[1, :], v[1, :], μ_v[1], μ_v[1])
    yy = get_covariance(v[2, :], v[2, :], μ_v[2], μ_v[2])
    zz = get_covariance(v[3, :], v[3, :], μ_v[3], μ_v[3])
    xy = get_covariance(v[2, :], v[1, :], μ_v[2], μ_v[1])
    xz = get_covariance(v[3, :], v[1, :], μ_v[3], μ_v[1])
    yz = get_covariance(v[3, :], v[2, :], μ_v[3], μ_v[2])
    return [
        xx xy xz;
        xy yy yz;
        xz yz zz
    ]
end
function get_covariance_mat(v::Matrix{Float64})
    xx = get_covariance(v[1, :], v[1, :])
    yy = get_covariance(v[2, :], v[2, :])
    zz = get_covariance(v[3, :], v[3, :])
    xy = get_covariance(v[2, :], v[1, :])
    xz = get_covariance(v[3, :], v[1, :])
    yz = get_covariance(v[3, :], v[2, :])
    return [
        xx xy xz;
        xy yy yz;
        xz yz zz
    ]
end


function get_linear_cd(x::Float64, y::Float64, μ::Float64, lin_fit)
    R2 = 1.0 - sum(abs2.(y - (lin_fit[1] * x .+ lin_fit[2])), dims = 1)[1] * sum(abs2.(y .- μ), dims = 1)[1]^-1
    return R2
end


function compensate_non_inertial(a::Vector{Float64}, φ::Float64, h::Float64)
    g = 9.780327 * (1.0 + 5.3024 * 10.0^-3 * sind(φ)^2 - 5.8 * 10.0^-6 * sind(2.0φ)^2) - (3.0877 * 10.0^-6 - 4.4 * 10.0^-9 * sind(φ)^2) * h + 7.2 * 10.0^-14 * h^2
    return a + [0.0; 0.0; g]
end


function get_eigenvalues_only(A::Matrix{Float64})
    return eigen(A).values
end


data = load(normpath(joinpath((@__FILE__), raw"..\AVT_data2.jld2")))


timestamp = data["times"]
timestamp .-= timestamp[1]
n_max = length(timestamp)


φ = data["lat"]
λ = data["lon"]
h = data["alt"]
v_x = data["v_x"]
v_y = data["v_y"]
v_z = data["v_z"]
GPS_valid = data["valid"]
φ′ = φ[GPS_valid]
λ′ = λ[GPS_valid]
h′ = h[GPS_valid]
v_x′ = v_x[GPS_valid]
v_y′ = v_y[GPS_valid]
v_z′ = v_z[GPS_valid]


ω_x = data["ω_x"]
ω_y = data["ω_y"]
ω_z = data["ω_z"]
a_x = data["a_x"]
a_y = data["a_y"]
a_z = data["a_z"]


m_x = data["m_x"]
m_y = data["m_y"]
m_z = data["m_z"]
HI = [40.953782; -14.121386; -36.173972]
SF = [
    1.378979 -0.036103 -0.069920;
    -0.036103 1.470051 -0.005797;
    -0.069920 -0.005797 1.455451
]
mag_valid = trues(n_max)
for i = 2:n_max
    if (m_x[i] == m_x[i - 1]) && (m_y[i] == m_y[i - 1]) && (m_z[i] == m_z[i - 1])
        mag_valid[i] = false
    end
end
m_x′ = m_x[mag_valid]
m_y′ = m_y[mag_valid]
m_z′ = m_z[mag_valid]
m_vec = [m_x′'; m_y′'; m_z′']
m_vec .-= HI
m_vec = SF * m_vec


p = data["pressure"]
h_p = p2alt.(p)
T = data["temp"] .+ 273.15


a = [a_x'; a_y'; a_z']
a_slices = [a[:, i:(i + win - 1)] for i = 1:(length(t) - win + 1)]
cov_a_slices = cov.(a_slices, dims = 2)
ev_a_slices = get_eigenvalues_only.(cov_a_slices)


#=
cov_GPS = get_covariance_mat([φ′'; λ′'; h′'])
cov_v = get_covariance_mat([v_x′'; v_y′'; v_z′'])


cov_ω = get_covariance_mat([ω_x'; ω_y'; ω_z'])
cov_a = get_covariance_mat([a_x'; a_y'; a_z'])


cov_m = get_covariance_mat([m_x′'; m_y′'; m_z′'])


μ_m = [
    get_mean(m_vec[1, :]);
    get_mean(m_vec[2, :]);
    get_mean(m_vec[3, :])
]
μ_m_normed = μ_m ./ norm(μ_m)
m_norms = sqrt.(sum(m_vec .* m_vec, dims = 1))
m_normed = m_vec ./ m_norms
mag_angular_dispersion = acos.(sum(m_normed .* μ_m_normed, dims = 1))
RMS_mag_angular_dispersion = sqrt(get_mean(abs2.(mag_angular_dispersion[:]))) * 180.0 / π


cov_p = get_covariance(p, p)
cov_h_p = get_covariance(h_p, h_p)
cov_p_T = get_covariance(p, T)
cov_T = get_covariance(T, T)


avg_GPS = [
    get_mean(φ);
    get_mean(λ);
    get_mean(h)
]


a_bias = [
    get_mean(a_x);
    get_mean(a_y);
    get_mean(a_z)
]
a_bias = compensate_non_inertial(a_bias, avg_GPS[1], avg_GPS[3])


ω_bias = [
    get_mean(ω_x);
    get_mean(ω_y);
    get_mean(ω_z)
]
=#




