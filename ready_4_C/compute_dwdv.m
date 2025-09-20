function dwdv = compute_dwdv(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude

[R_phi, R_lamb, ~, ~] = compute_radii(x);

% Eqn 7.74b
m12 = 1/(R_lamb + h);
m21 = -1/(R_phi + h);
m32 = -tand(phi)/(R_lamb + h);

% assemble dwdv matrix
dwdv = [0, m12, 0;...
    m21, 0, 0;...
    0, m32, 0];

end