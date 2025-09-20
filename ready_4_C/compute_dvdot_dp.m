function dvdot_dp = compute_dvdot_dp(x,we)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); vd = x(10); % get vel components

[R_phi, R_lamb, dR_phi_dphi, dR_lamb_dphi] = compute_radii(x);

[~,dg_dphi,dg_dh] = compute_g_dg(x);

Y11 = -(ve^2*secd(phi)^2)/(R_lamb + h) + (ve^2*tand(phi))/(R_lamb + h)^2*dR_lamb_dphi ... % Eqn 7.82a
    - 2*we*ve*cosd(phi) - (vn*vd)/(R_phi + h)^2 * dR_phi_dphi;
Y13 = (ve^2*tand(phi))/(R_lamb + h)^2 - (vn*vd)/(R_phi + h)^2; % Eqn 7.82b
Y21 = (ve*vn*secd(phi)^2)/(R_lamb + h) - (ve*vn*tand(phi))/(R_lamb + h)^2*dR_lamb_dphi + ... % Eqn 7.82c
    2*we*vn*cosd(phi) - (ve*vd)/(R_lamb + h)^2*dR_lamb_dphi - 2*we*vd*sind(phi);
Y23 = -ve*((vn*tand(phi) + vd)/(R_lamb + h)^2); % Eqn 7.82d
Y31 = ve^2/(R_lamb + h)^2*dR_lamb_dphi + vn^2/(R_phi + h)^2*dR_phi_dphi + 2*we*ve*sind(phi) + dg_dphi; % Eqn 7.82e
Y33 = ve^2/(R_lamb + h)^2 + vn^2/(R_phi + h)^2 + dg_dh; % Eqn 7.82f

% assemble dvdot_dp matrix
dvdot_dp = [Y11, 0, Y13;...
    Y21, 0, Y23;...
    Y31, 0, Y33];
end