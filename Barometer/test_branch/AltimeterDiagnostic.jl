using LinearAlgebra
using Statistics


using SatelliteToolboxAtmosphericModels
using JLD2
using GLMakie
using BenchmarkTools


include("AltimeterPressureTesting.jl")


begin
    # FAR coordinates: 35.34788, -117.80683, 625.0
    t_launch = DateTime("2025-01-10T11:00:0.0")
    φ_launch = 35.34788
    λ_launch = -117.80683
    h_launch = 625.0

    # molar masses from PubChem
    local molar_masses::Vector{Float64} = [
        0.01400674;
        0.028014;
        0.0159994;
        0.031999;
        0.00100794;
        0.004002602; 
        0.039948
    ]
    

    local function normalized_abundances(atm_data)
        local Σ_number_density = (
            atm_data.N_number_density +
            atm_data.N2_number_density +
            atm_data.O_number_density + atm_data.aO_number_density +
            atm_data.O2_number_density +
            atm_data.H_number_density +
            atm_data.He_number_density + 
            atm_data.Ar_number_density
        )
        return abundances = Σ_number_density^-1 * [
            atm_data.N_number_density;
            atm_data.N2_number_density;
            atm_data.O_number_density + atm_data.aO_number_density;
            atm_data.O2_number_density;
            atm_data.H_number_density;
            atm_data.He_number_density;
            atm_data.Ar_number_density
        ]
    end


    """
        pressure_nrlmsise00(h)

    Return the barometric pressure of the atmosphere at the specified geodetic altitude.

    ...
    # Arguments
    - `h::Real`: Geodetic altitude (metres)

    ...
    """
    function pressure_nrlmsise00(h)
        local atm_data = AtmosphericModels.nrlmsise00(
            t_launch,
            h,
            φ_launch,
            λ_launch
        )
        local frac_species = normalized_abundances(atm_data)
        local T = atm_data.temperature
        local m = dot(molar_masses, frac_species)
        local ρ = atm_data.total_density
        return p = ρ * (8.31446261815324 * m^-1 * T)
    end


    """
        true_log_pressure(h)

    Return the natural logarithm of the barometric pressure of the atmosphere at the specified geodetic altitude.

    ...
    # Arguments
    - `h::Real`: Geodetic altitude (metres)

    ...
    """
    function true_log_pressure(h)
        return log(pressure_nrlmsise00(h))
    end


    """
        true_log_deriv_pressure(h)

    Return the logarithmic derivative with respect to geodetic altitude of the barometric pressure of the atmosphere at the specified geodetic altitude.

    ...
    # Arguments
    - `h::Real`: Geodetic altitude (metres)

    ...
    """
    function true_log_deriv_pressure(h)
        return ForwardDiff.derivative(true_log_pressure, h)
    end


    """
        true_log_deriv_pressure(h)

    Return the second derivative with respect to geodetic altitude of the logarithm of the barometric pressure of the atmosphere at the specified geodetic altitude.

    ...
    # Arguments
    - `h::Real`: Geodetic altitude (metres)

    ...
    """
    function true_log_2nd_deriv_pressure(h)
        return ForwardDiff.derivative(true_log_deriv_pressure, h)
    end
end


h_min::BigFloat = -1.0 * 10.0^3
h_max::BigFloat = 12.0 * 10.0^3
h_samples::Vector{BigFloat} = collect(range(h_min, h_max, length = 2^12))
true_pressure::Vector{BigFloat} = pressure_nrlmsise00.(h_samples)


h_1::Vector{Float32} = pressure_altimeter_uncorrected.(Float32.(true_pressure), (1, ))
h_2::Vector{Float32} = pressure_altimeter_uncorrected.(Float32.(true_pressure), (2, ))
h_4::Vector{Float32} = pressure_altimeter_uncorrected.(Float32.(true_pressure), (4, ))
h_6::Vector{Float32} = pressure_altimeter_uncorrected.(Float32.(true_pressure), (6, ))
h_8::Vector{Float32} = pressure_altimeter_uncorrected.(Float32.(true_pressure), (8, ))
h_12::Vector{Float32} = pressure_altimeter_uncorrected.(Float32.(true_pressure), (12, ))
h_16::Vector{Float32} = pressure_altimeter_uncorrected.(Float32.(true_pressure), (16, ))


ε_1::Vector{BigFloat} = abs.(BigFloat.(h_1) - h_samples)
ε_2::Vector{BigFloat} = abs.(BigFloat.(h_2) - h_samples)
ε_4::Vector{BigFloat} = abs.(BigFloat.(h_4) - h_samples)
ε_6::Vector{BigFloat} = abs.(BigFloat.(h_6) - h_samples)
ε_8::Vector{BigFloat} = abs.(BigFloat.(h_8) - h_samples)
ε_12::Vector{BigFloat} = abs.(BigFloat.(h_12) - h_samples)
ε_16::Vector{BigFloat} = abs.(BigFloat.(h_16) - h_samples)


mean_errors = [
    sqrt(mean(abs2.(ε_1)));
    sqrt(mean(abs2.(ε_2)));
    sqrt(mean(abs2.(ε_4)));
    sqrt(mean(abs2.(ε_6)));
    sqrt(mean(abs2.(ε_8)));
    sqrt(mean(abs2.(ε_12)));
    sqrt(mean(abs2.(ε_16)))
]
max_errors = [
    maximum(ε_1);
    maximum(ε_2);
    maximum(ε_4);
    maximum(ε_6);
    maximum(ε_8);
    maximum(ε_12);
    maximum(ε_16)
]


fig = Figure(size = (4000, 2000), px_per_unit = 1)
ax = Axis(fig[1, 1], title = "Laguerre method truncation error", xlabel = "true altitude (metres)", ylabel = "absolute error (metres)")
lines!(ax, true_pressure, ε_1, label = "1 iteration")
lines!(ax, true_pressure, ε_2, label = "2 iterations")
lines!(ax, true_pressure, ε_16, label = "16 iterations")
axislegend(ax, position = :lt)
fig
save(normpath(joinpath((@__FILE__), raw"..\altimeter_test_log_plot.png")), fig)
#display(fig, px_per_unit = 2)




#typeof(trial_output)

