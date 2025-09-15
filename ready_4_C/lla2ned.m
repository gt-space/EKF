function ned = lla2ned(lla,lla0)
ecef = lla2ecef(lla);
ned = ecef2ned(ecef,lla0);
end