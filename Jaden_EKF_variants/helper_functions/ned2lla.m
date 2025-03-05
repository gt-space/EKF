function lla = ned2lla(ned,lla0)

lat0 = deg2rad(lla0(1));
lon0 = deg2rad(lla0(2));
h0 = lla0(3);

% WGS84 parameters
a = 6378137;           % semi-major axis (m)
e = 0.08181919;        % eccentricity

% Radii of curvature at the origin
N_val = a / sqrt(1 - e^2 * sin(lat0)^2);
M_val = a*(1 - e^2) / (1 - e^2 * sin(lat0)^2)^(3/2);

% Given NED offset (in meters)
n = ned(1);   % e.g., 100 m north
e_offset = ned(2);  % e.g., 50 m east
d = ned(3);    % e.g., 10 m down

% Compute new LLA (in radians for lat, lon)
lat_new = lat0 + n / (M_val + h0);
lon_new = lon0 + e_offset / ((N_val + h0) * cos(lat0));
h_new = h0 - d;

% Convert lat and lon back to degrees
lat_new_deg = rad2deg(lat_new);
lon_new_deg = rad2deg(lon_new);
lla = [lat_new_deg;lon_new_deg;h_new];
end