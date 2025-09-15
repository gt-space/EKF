function ecef = ned2ecef(ned,lla0)

ecef0 = lla2ecef(lla0); % get ECEF of NED origin
phi = lla0(1);
lambda = lla0(2);
s1 = sind(phi);
c1 = cosd(phi);
s2 = sind(lambda);
c2 = cosd(lambda);
% compute NED2ECEF DCM
An2e = [-s1*c2 -s2 -c1*c2;
        -s1*s2 c2  -c1*s2;
         c1    0      -s1];
ecef = double(ecef0 + An2e*ned);
end