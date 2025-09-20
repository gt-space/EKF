function dpdot_dp = compute_dpdot_dp(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); % get vel components

[R_phi, R_lamb, dR_phi_dphi, dR_lamb_dphi] = compute_radii(x);

% Eqn 7.80a
m11 = -vn/(R_phi + h)^2*dR_phi_dphi;
m13 = -vn/(R_phi + h)^2;
m21 = -(ve*secd(phi))/(R_lamb + h)^2*dR_lamb_dphi + (ve*secd(phi)*tand(phi))/(R_lamb + h);
m23 = -ve*secd(phi)/(R_lamb + h)^2;

% assemble dpdot_dp matrix
dpdot_dp = [m11, 0, m13;...
    m21, 0, m23;...
    0, 0, 0];

end