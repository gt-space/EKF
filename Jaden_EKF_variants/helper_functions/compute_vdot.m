function vdot = compute_vdot(x,ahat_n)
we = 7.2921158e-5; % rad/s, Earth rotation rate
phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); vd = x(10); % get vel components
an = ahat_n(1); ae = ahat_n(2); ad = ahat_n(3); % get accel components

[R_phi, R_lamb, ~, ~] = compute_radii(x); % compute radii
[ghat,~,~] = compute_g(x); % get gravity estimate

vndot = -((ve/(R_lamb + h)*cosd(phi)) + 2*we)*ve*sind(phi) + (vn*vd)/(R_phi + h) + an; % Eqn 7.68f
vedot = ((ve/(R_lamb + h)*cosd(phi)) + 2*we)*vn*sind(phi) + (ve*vd)/(R_lamb + h) + 2*we*vd*cosd(phi) + ae; % Eqn 7.68g
vddot = -ve^2/(R_lamb + h) - vn^2/(R_phi + h) - 2*we*ve*cosd(phi) + ghat + ad; % Eqn 7.68h

vdot = [vndot; vedot; vddot]; % assemble vdot vector
end