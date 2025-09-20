function wn = compute_wn(x,we)
% see Eq A.272 for why we need w N to I expressed in N frame coordinates
% N is NED, I is Earth-Center Inertial, see A.9.1

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); % north velocity
ve = x(9); % east velocity

[R_phi, R_lamb, ~, ~] = compute_radii(x);

vec1 = we*[cosd(phi); 0; -sind(phi)];
vec2 = [ve/(R_lamb + h); -vn/(R_phi + h); -(ve*tand(phi))/(R_lamb + h)];

wn = vec1 + vec2; % Eqn 7.70

end