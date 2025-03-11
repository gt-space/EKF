function lla = ecef2lla(ecef)
a = 6378137;
b = 6356752.3124;
e = 0.08181919;
x = ecef(1);
y = ecef(2);
z = ecef(3);
p = norm(ecef(1:2));
psi = atan2(z*a,(p*b));
ebar2 = (a/b)^2-1;
phi = rad2deg(atan2(z+ebar2*b*sin(psi)^3,(p-e^2*a*cos(psi)^3)));
N = a/sqrt(1-e^2*sind(phi)^2);
lambda = rad2deg(atan2(y,x));
h = p/cosd(phi) - N;
lla = [phi;lambda;h];
end