% derivative of the pressure function
function dpdh = pressure_derivative(h, p_0, a)
    dpdh = p_0 * (a * [0.0; 1.0; 2.0 * h; 3.0 * h^2; 4.0 * h^3; 5.0 * h^4; 6.0 * h^5; 7.0 * h^6; 8.0 * h^7; 9.0 * h^8; 10.0 * h^9; 11.0 * h^10; 12.0 * h^11; 13.0 * h^12; 14.0 * h^13; 15.0 * h^14]) * 10^(a * [1.0; h; h^2; h^3; h^4; h^5; h^6; h^7; h^8; h^9; h^10; h^11; h^12; h^13; h^14; h^15]);
    dpdh = dpdh * 10^-5
end
