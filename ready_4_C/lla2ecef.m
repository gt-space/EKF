function ecef = lla2ecef(lla)
% See book ""Optimal Estimation of Dynamic Systems" Appendix A
% verified with Google Maps
a = 6378137; % semi-major axis
b = 6356752.3124; % semi-minor axis
e = 0.08181919; % eccentricity
phi = lla(1);
lambda = lla(2); 
h = lla(3);
N = a/sqrt(1-(e*sind(phi))^2); % A.238
x = (N+h)*cosd(phi)*cosd(lambda); % A.239
y = (N+h)*cosd(phi)*sind(lambda);
z = (N*(1-e^2)+h)*sind(phi);
ecef = double([x;y;z]);
end