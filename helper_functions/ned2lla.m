function lla = ned2lla(ned,lla0)
ecef = ned2ecef(ned,lla0);
lla = ecef2lla(ecef);
end