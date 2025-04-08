load('raw_flight_sim_data.mat');
wgs84 = wgs84Ellipsoid('meter');

posI = simData.pos;
lla0 = simData.lla0;
lat0 = lla0(1);
lon0 = lla0(2);
h0 = lla0(3);
% n = posI(1,1);
% e = posI(2,1);
% d = posI(3,1);
% 
% lla0_old = ned2lla(posI(:,1),lla0);
% [lat_true,lon_true,h_true] = ned2geodetic(n,e,d,lat0,lon0,h0,wgs84,"degrees");

% ecef = ned2ecef(posI(:,1),lla0);
% lla0_new = improved_ecef2lla(ecef);

% [X,Y,Z] = geodetic2ecef(wgs84,lat,lon,h);
% 
% ecef_true = [X;Y;Z];

lla_old = zeros(size(posI));
lla_true = zeros(size(posI));

% lla_new = zeros(size(posI));
for i = 1:size(posI,2)
    lla_old(:,i) = ned2lla(posI(:,i),lla0);
    [lat_true,lon_true,h_true] = ned2geodetic(posI(1,i),posI(2,i),posI(3,i),lat0,lon0,h0,wgs84,"degrees");
    lla_true(:,i) = [lat_true;lon_true;h_true];
    ecef = ned2ecef(posI(:,i),lla0);
    lla_new(:,i) = improved_ecef2lla(ecef);
end
figure
plot((lla_old(1,1:100:end)-lla_true(1,1:100:end)))
figure
plot((lla_old(2,1:100:end)-lla_true(2,1:100:end)))
figure
plot(lla_old(3,1:100:end)-lla_true(3,1:100:end))