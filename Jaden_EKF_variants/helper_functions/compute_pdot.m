function pdot = compute_pdot(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); % north velocity
ve = x(9); % east velocity
vd = x(10); % down velocity

[R_phi, R_lamb, ~, ~] = compute_radii(x); % compute radii

phidot = vn/(R_phi + h); % Eqn 7.68c
lambdot = ve/((R_lamb + h)*cosd(phi)); % Eqn 7.68d
hdot = -vd; % Eqn 7.68e

pdot = [phidot; lambdot; hdot]; % assemble pdot vector
end