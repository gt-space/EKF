clear;
load("raw_flight_sim_data.mat");
% apogee_index = 9447;
final_index = 10000; % 100 secs, covering apogee
posI = simData.pos(:,1:final_index); % NED
accI = simData.acc(:,1:final_index); % NED
q = double(simData.quat(:,1:final_index));
omegaB = simData.omega(:,1:final_index); % body
time = double(simData.time(:,1:final_index));
dt = time(2) - time(1);
IMU_freq = round(1/dt); % Hz
accB = zeros(size(accI));
for i = 1:length(time)
    accB(:,i) = q2v(sandwich(qconj(q(:,i)),[0;accI(:,i)])); % I2B
end
%% generate IMU outputs
acc = zeros(size(accB)); % IMU output of acceleration
nla = 5e-3; % nonlinearity of accelerometer, 0.5% of +-16g  
ba = zeros(size(omegaB)); % bias of accelerometer
Sa = ones(size(omegaB)); % scaling factor of accelerometer
gyro = zeros(size(omegaB)); % IMU output of angular velocity
nlg = 1e-4; % nonlinearity of gyro, 0.01% of +-2000 dps
bg = zeros(size(omegaB)); % bias of gyro
Sg = ones(size(omegaB)); % scaling factor of gyro
% random walk TODO
for i = 2:length(time)
    Sa(:,i) = Sa(:,i-1) + nla/10*(rand(3,1)-0.5);
    ba(:,i) = ba(:,i-1) + dt/10*40e-6 *9.8*randn(3,1);   % 40 mug
    Sg(:,i) = Sg(:,i-1) + nlg/10*(rand(3,1)-0.5);
    bg(:,i) = bg(:,i-1) + dt/10*deg2rad(3/3600)*randn(3,1);   
end
R_acc = eye(3) * (200e-6*9.81)*sqrt(IMU_freq); % std of acc noise, 200 mug/sqrt(Hz)
R_gyro = eye(3) * deg2rad(12e-3)*sqrt(IMU_freq); % std of gyro noise, 12 mdps/sqrt(Hz)
acc_noise = R_acc*randn(size(acc));
gyro_noise = R_gyro*randn(size(gyro));

for i = 1:length(time)
    acc(:,i) = Sa(:,i) .* accB(:,i) + ba(:,i) + acc_noise(:,i); 
    gyro(:,i) = Sg(:,i) .* omegaB(:,i) + bg(:,i) + gyro_noise(:,i); 
end

%% generate GPS data
lla0 = simData.lla0;
GPS_freq = 2; % Hz
ratio = round(IMU_freq / GPS_freq);
LLA = zeros(3, length(time));
for i = 1: length(time)
    LLA(:,i) = ned2lla(posI(:,i),lla0); % need to verify ned2lla TODO
end
R_gps = diag([1.35e-5;1.63e-5;2]); % GPS noise TODO
noisy_LLA = (LLA + R_gps * randn(size(LLA)));
GPS = nan(size(LLA));
for i = 1: ratio: length(time)
    GPS(:,i) = noisy_LLA(:,i); 
end
% GPS(:,1:3000) = nan(3,3000); % simulated GPS cutoff
% GPS(:,6001:end) = nan(3,4000);
%% generate magnotometer data
% TODO need sensor datasheet
magI = [0.4891; 0.104; 0.866]; % need more precise vector at launch site TODO
mag_freq = 10; % Hz
ratio = round(IMU_freq / mag_freq);
R_mag = diag([3.2e-7 4.1e-7 3.2e-7]); % magnetometer noise TODO
magB = nan(size(bg));
mag_noise = R_mag * randn(size(magB));
for i = 1:ratio:size(magB,2)
    magB(:,i) = quaternion2DCM(q(:,i))' * magI ;
end
magB = magB + mag_noise;
%%
sim_data.time = time;
sim_data.acc = acc;
sim_data.gyro = gyro; % TODO what is the output of gyro? Euler rate?
sim_data.ba = ba; % "true" bias
sim_data.bg = bg;
sim_data.Sa = Sa; % "true" scaling factor
sim_data.Sg = Sg;
sim_data.LLA = LLA; % "true" LLA position
sim_data.q = q; % "true" attitude
sim_data.GPS = GPS; % GPS measurment
sim_data.magB = magB; % magnetometer measurement
sim_data.magI = magI;
save("sim_data_4_EKF_testing.mat","sim_data");
