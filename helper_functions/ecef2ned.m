function ned = ecef2ned(ecef,lla0)

ecef0 = lla2ecef(lla0); % get ECEF of NED origin
phi = lla0(1);
lambda = lla0(2);
s1 = sind(phi);
c1 = cosd(phi);
s2 = sind(lambda);
c2 = cosd(lambda);
% get ECEF2NED DCM, A.242
Ae2n = [-s1*c2 -s1*s2  c1;
        -s2     c2     0;
        -c1*c2  -c1*s2 -s1];
ned = double(Ae2n*(ecef - ecef0));
end