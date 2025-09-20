function dvdot_dv = compute_dvdot_dv(x,we)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); vd = x(10); % get vel components

[R_phi, R_lamb,~,~] = compute_radii(x);

Z11 = vd/(R_phi + h); % Eqn 7.83a
Z12 = (-2*ve*tand(phi))/(R_lamb + h) + 2*we*sind(phi); % Eqn 7.83a
Z13 = vn/(R_phi + h); % Eqn 7.83a
Z21 = (ve*tand(phi))/(R_lamb + h) + 2*we*sind(phi); % Eqn 7.83b
Z22 = (vd + vn*tand(phi))/(R_lamb + h); % Eqn 7.83b
Z23 = ve/(R_lamb + h) + 2*we*cosd(phi); % Eqn 7.83b
Z31 = (-2*vn)/(R_phi + h); % Eqn 7.83c
Z32 = (-2*ve)/(R_lamb + h) - 2*we*cosd(phi); % Eqn 7.83c

% assemble dvdot_dv matrix
dvdot_dv = [Z11, Z12, Z13;...
    Z21, Z22, Z23;...
    Z31, Z32,  0];

end