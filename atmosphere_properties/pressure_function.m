% uses 7th order polynomial best fit
% height in meters
% pressure in pascals
% turns an altitude, h (meters), into a pressure, p
function p = pressure_function(h, p_0, a)
    p = 10^(5.0122185 + (-4.9929004 * 10^-5 * h) + (-5.415637 * 10^-10 * h^2) + (-3.837231 * 10^-14 * h^3) + (2.55155 * 10^-18 * h^4) + (-5.321706 * 10^-23 * h^5) + (4.813401 * 10^-28 * h^6) + (-1.6294356 * 10^-33 * h^7))
end


% 0.9999895 R^2 from -2 km to 80 km
% simulated for 35.3478813171387, -117.806831359863 (lat, long) at 2025-04-25 14:00:00 UTC 



