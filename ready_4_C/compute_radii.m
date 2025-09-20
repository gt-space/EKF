function [R_phi, R_lamb, dR_phi_dphi, dR_lamb_dphi] = compute_radii(x)

% WGS84 model parameters (source: https://en.wikipedia.org/wiki/World_Geodetic_System)
a = 6378137; % semi-major axis [m]
b = 6356752.31425; % semi-minor axis [m]
e2 = 1 - (b/a)^2; % eccentricity

phi = x(5); % get latitude

num1 = a*(1-e2);
den1 = (1-e2*sind(phi)^2)^(3/2);

den2 = (1-e2*sind(phi)^2)^(1/2);

R_phi = num1/den1; % Eqn 7.69a
R_lamb = a/den2; % Eqn 7.69b

num3 = 3*a*(1-e2)*e2*sind(phi)*cosd(phi);
den3 = (1-e2*sind(phi)^2)^(5/2);

num4 = a*e2*sind(phi)*cosd(phi);
den4 = (1-e2*sind(phi)^2)^(3/2);

dR_lamb_dphi = num4/den4; % Eqn 7.75a
dR_phi_dphi = num3/den3; % Eqn 7.75b

end