% derivative of the pressure function
function dpdh = pressure_derivative(h, p_0, a)
    dpdh = log(10.0) * (a(2) + 2.0 * a(3) * h + 3.0 * a(4) * h^2 + 4.0 * a(5) * h^3 + 5.0 * a(6) * h^4 + 6.0 * a(7) * h^5 + 7.0 * a(8) * h^6) * 10^(a(1) + a(2) * h + a(3) * h^2 + a(4) * h^3 + a(5) * h^4 + a(6) * h^5 + a(7) * h^6 + a(8) * h^7)
end
