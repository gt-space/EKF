function dwdp = compute_dwdp(x,we)
%Eqn 7.74a

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); % get vel components

[R_phi, R_lamb, dR_phi_dphi, dR_lamb_dphi] = compute_radii(x);

% Eqn 7.74a
m11 = -we*sind(phi) - ve/(R_lamb + h)^2 * dR_lamb_dphi;
m13 = -ve/(R_lamb + h)^2;
m21 = vn/(R_phi + h)^2 * dR_phi_dphi;
m23 = vn/(R_phi + h)^2;
m31 = -we*cosd(phi)-(ve*secd(phi)^2)/(R_lamb + h) + (ve*tand(phi))/(R_lamb + h)^2*dR_lamb_dphi;
m33 = (ve*tand(phi))/(R_lamb + h)^2;

% assemble dwdp matrix
dwdp = [m11, 0, m13;...
    m21, 0, m23;...
    m31, 0, m33];

end