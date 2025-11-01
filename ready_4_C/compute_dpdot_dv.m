function dpdot_dv = compute_dpdot_dv(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude

[R_phi, R_lamb, ~, ~] = compute_radii(x);

% Eqn 7.80b
m11 = rad2deg(  1/(R_phi + h) );
m22 = rad2deg( secd(phi)/(R_lamb + h) );
m33 = -1;

% assemble dpdot_dv matrix
dpdot_dv = diag([m11, m22, m33]);
end