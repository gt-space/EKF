% uses 7th order polynomial best fit (first entry in a is constant term, then second entry is 1st order, third entry is 2nd order, etc.)
% height in meters
% pressure in pascals
% turns an altitude, h (meters), into a pressure, p
function p = pressure_function(h, p_0, a)
    p = 10^(a(1) + a(2) * h + a(3) * h^2 + a(4) * h^3 + a(5) * h^4 + a(6) * h^5 + a(7) * h^6 + a(8) * h^7)
end


% 0.9999895 R^2 from -2 km to 80 km
% simulated for 35.3478813171387, -117.806831359863 (lat, long) at 2025-04-25 14:00:00 UTC 
%5.0122185
%-4.9929004e-5
%-5.415637e-10
%-3.837231e-14
%2.55155e-18
%-5.321706e-23
%4.813401e-28
%-1.6294356e-33




