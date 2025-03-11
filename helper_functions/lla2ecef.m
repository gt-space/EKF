function ecef = lla2ecef(lla)
a = 6378137;
b = 6356752.3124;
e = 0.08181919;
phi = lla(1);
lambda = lla(2);
h = lla(3);
N = a/sqrt(1-(e*sind(phi))^2);
x = (N+h)*cosd(phi)*cosd(lambda);
y = (N+h)*cosd(phi)*sind(lambda);
z = (N*(1-e^2)+h)*sind(phi);
ecef = double([x;y;z]);
end