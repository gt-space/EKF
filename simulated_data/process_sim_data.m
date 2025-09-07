clear;
load("raw_flight_sim_data.mat");
% apogee_index = 9447;
apogee_index = find(simData.vel(3,:)>0,1);
final_index = length(simData.time); % 100 secs, covering apogee
% final_index = 104932; % 100 secs, covering apogee
% final_index = apogee_index+30000; % 100 secs, covering apogee
posI = double(simData.pos(:,1:final_index)); % NED
accI = double(simData.acc(:,1:final_index)); % NED
q = double(simData.quat(:,1:final_index));
omegaB = double(simData.omega(:,1:final_index)); % body
time = double(simData.time(:,1:final_index));
dt = time(2) - time(1);
IMU_freq = round(1/dt); % Hz
%% generate GPS data
lla0 = double(simData.lla0);
GPS_freq = 1; % Hz
ratio = round(IMU_freq / GPS_freq);
LLA = zeros(3, length(time));
for i = 1: length(time)
    LLA(:,i) = ned2lla(posI(:,i),lla0); 
end
R_gps = diag([1.35e-5;1.65e-5;2]); % GPS noise TODO
noisy_LLA = (LLA + R_gps * randn(size(LLA)));
GPS = nan(size(LLA));
for i = 1: ratio: length(time)
    if noisy_LLA(3,i) < 29000 && norm(accI(:,i)) < 4*9.8
        GPS(:,i) = noisy_LLA(:,i);
    end
end
% GPS(:,1:3000) = nan(3,3000); % simulated GPS cutoff
% GPS(:,6001:end) = nan(3,4000);
%%
accB = zeros(size(accI));
for i = 1:length(time)
    ghat = [0;0;compute_g(LLA(:,i))];
    accB(:,i) = q2v(sandwich(qconj(q(:,i)),[0;accI(:,i)-ghat])); % I2B
end
%% generate IMU outputs
acc = zeros(size(accB)); % IMU output of acceleration
nla = 5e-4; % nonlinearity of accelerometer, 0.5% of +-16g  
ba = zeros(size(omegaB)); % bias of accelerometer
Sa = ones(size(omegaB)); % scaling factor of accelerometer
gyro = zeros(size(omegaB)); % IMU output of angular velocity
nlg = 1e-5; % nonlinearity of gyro, 0.01% of +-2000 dps
bg = zeros(size(omegaB)); % bias of gyro
Sg = ones(size(omegaB)); % scaling factor of gyro
sigv_a = 200e-6 * 9.81;  % std of acc noise, 200 mug/sqrt(Hz)
sigv_g = deg2rad(12e-3); % std of gyro noise, 12 mdps/sqrt(Hz)
sigu_a = 40e-6 * 9.81; % acc bias instability, 40 mug
sigu_g = deg2rad(3/3600); % gyro bias instability, 3 degree per hour
% random walk TODO
for i = 2:length(time)
    Sa(:,i) = Sa(:,i-1) + nla/10*(rand(3,1)-0.5);
    ba(:,i) = ba(:,i-1) + sqrt(dt)*sigu_a*randn(3,1);   % 40 mug
    Sg(:,i) = Sg(:,i-1) + nlg/10*(rand(3,1)-0.5);
    bg(:,i) = bg(:,i-1) + sqrt(dt)*sigu_g*randn(3,1);   
end

acc(:,1) = Sa(:,1) .* accB(:,1) + 1/2*ba(:,1) + sqrt(sigv_a^2/dt + 1/12*sigu_a^2 * dt)* randn(3,1); 
gyro(:,1) = Sg(:,1) .* omegaB(:,1) + 1/2*bg(:,1) + sqrt(sigv_g^2/dt + 1/12*sigu_g^2 * dt)* randn(3,1);
for i = 2:length(time)
    acc(:,i) = Sa(:,i) .* accB(:,i) + 1/2*(ba(:,i)+ba(:,i-1)) + sqrt(sigv_a^2/dt + 1/12*sigu_a^2 * dt)* randn(3,1); 
    gyro(:,i) = Sg(:,i) .* omegaB(:,i) + 1/2*(bg(:,i)+bg(:,i-1)) + sqrt(sigv_g^2/dt + 1/12*sigu_g^2 * dt)* randn(3,1);
end


%% generate magnotometer data
% TODO need sensor datasheet
magI = [0.4891; 0.104; 0.866]; % need more precise vector at launch site TODO
mag_freq = 10; % Hz
ratio = round(IMU_freq / mag_freq);
R_mag = diag([3.2e-7 4.1e-7 3.2e-7]); % magnetometer noise TODO in Tesla
magB = nan(size(bg));
mag_noise = R_mag * randn(size(magB));
for i = 1:ratio:size(magB,2)
    magB(:,i) = quaternion2DCM(q(:,i))' * magI ;
end
magB = magB + mag_noise;
%% Generate Barometer Measurements
baro_freq = 20;
ratio = round(IMU_freq / baro_freq);
p_0 = 102715.47296217596;
a = [-0.004366921037228173;
   -5.0964993627948606 * 10^-5;
   1.6319216411827967 * 10^-9;
   -6.053135304413484 * 10^-13;
   6.911049045104908 * 10^-17;
   -4.516216131249696 * 10^-21;
   1.90623861182507 * 10^-25;
   -5.457995055935324 * 10^-30;
   1.0906415917824724 * 10^-34;
   -1.5456197431929747 * 10^-39;
   1.5609417634668778 * 10^-44;
   -1.1149608634620734 * 10^-49;
   5.500087730582754 * 10^-55;
   -1.7817554862391996 * 10^-60;
   3.409963106133886 * 10^-66;
   -2.9213421378235885 * 10^-72];
press = nan(size(time));
for i = (apogee_index + 5000) : ratio: final_index
    h = LLA(3,i);
    press(i) = pressure_function(h,p_0,a);
end
sig_baro =2.5e-3; % bar
baro_noise = sig_baro*randn(size(press));
baro = press + baro_noise;
%%
sim_data.time = time;
sim_data.acc = acc;
sim_data.gyro = gyro; % 
sim_data.ba = ba; % "true" bias
sim_data.bg = bg;
sim_data.Sa = Sa; % "true" scaling factor
sim_data.Sg = Sg;
sim_data.LLA0 = lla0; % "true" NED position
sim_data.NED = posI; % "true" NED position
sim_data.velI = double(simData.vel);
sim_data.q = q; % "true" attitude
sim_data.GPS = GPS; % GPS measurment
sim_data.magB = magB; % magnetometer measurement
sim_data.magI = magI;
sim_data.apogee_idx = apogee_index;
sim_data.baro = baro;
sim_data.press = press;
save("sim_data_4_EKF_testing.mat","sim_data");

%%
function [g,dg_dphi,dg_dh] = compute_g(x)

phi = x(1); % get latitude
h = x(3); % get height/altitude
g = 9.780327*(1 + 5.3024*10^(-3)*sind(phi)^2 - 5.8*10^(-6)*sind(2*phi)^2) - ...
     (3.0877*10^(-6) - 4.4*10^(-9)*sind(phi)^2)*h + 7.2*10^(-14)*h^2; % Eqn 7.69c
% g = 9.81;

dg_dphi = 9.780327*(1.06048*10^(-2)*sind(phi)*cosd(phi)... % Eqn 7.84a
    - 4.64*10^(-5)*(sind(phi)*cosd(phi)^3 - sind(phi)^3*cosd(phi)))...
    + 8.8*10^(-9)*h*sind(phi)*cosd(phi);

dg_dh = -3.0877*10^(-6) + 4.4*10^(-9)*sind(phi)^2 + 1.44*10^(-13)*h; % Eqn 7.84b

end