function [g,dg_dphi,dg_dh] = compute_g(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude
gWGS0 = 9.78032677;
g = gWGS0*(1 + 5.3024*10^(-3)*sind(phi)^2 - 5.8*10^(-6)*sind(2*phi)^2) - ...
     (3.0877*10^(-6) - 4.4*10^(-9)*sind(phi)^2)*h + 7.2*10^(-14)*h^2; % Eqn 7.69c
% g = 9.81;

dg_dphi = gWGS0*(1.06048*10^(-2)*sind(phi)*cosd(phi)... % Eqn 7.84a
    - 4.64*10^(-5)*(sind(phi)*cosd(phi)^3 - sind(phi)^3*cosd(phi)))...
    + 8.8*10^(-9)*h*sind(phi)*cosd(phi);

dg_dh = -3.0877*10^(-6) + 4.4*10^(-9)*sind(phi)^2 + 1.44*10^(-13)*h; % Eqn 7.84b

end