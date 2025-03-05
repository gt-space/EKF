clear;
load("MCrun.mat");
apogee_index = 9447;
posI = simData.pos(:,1:10000);
accI = simData.acc(:,1:10000);
q = double(simData.quat(:,1:10000));
omegaB = simData.omega(:,1:10000);
time = double(simData.time(:,1:10000));
dt = time(2) - time(1);
IMU_freq = round(1/dt); % Hz
accB = zeros(size(accI));
for i = 1:length(time)
    accB(:,i) = q2v(sandwich(qconj(q(:,i)),[0;accI(:,i)]));
end

acc = zeros(size(accB));
nla = 5e-3; % 0.5% of +-16g  
ba = zeros(size(omegaB));
Sa = ones(size(omegaB));
gyro = zeros(size(omegaB));
nlg = 1e-4; % 0.01% of +-2000 dps
bg = zeros(size(omegaB));
Sg = ones(size(omegaB));
for i = 2:length(time)
    Sa(:,i) = Sa(:,i-1) + nla*2*(rand(3,1)-0.5);
    ba(:,i) = ba(:,i-1) + dt*40e-6 *9.8*randn(3,1);   % 40 mug
    Sg(:,i) = Sg(:,i-1) + nlg*2*(rand(3,1)-0.5);
    bg(:,i) = bg(:,i-1) + dt*deg2rad(3/3600)*randn(3,1);   
end
S_acc = eye(3) * (200e-6*9.81)*sqrt(IMU_freq); % 200 mug/sqrt(Hz)
S_gyro = eye(3) * deg2rad(12e-3)*sqrt(IMU_freq); % 12 mdps/sqrt(Hz)
acc_noise = S_acc*randn(size(acc));
gyro_noise = S_gyro*randn(size(gyro));

for i = 1:length(time)
    acc(:,i) = Sa(:,i) .* accB(:,i) + ba(:,i) + acc_noise(:,i); 
    gyro(:,i) = Sg(:,i) .* omegaB(:,i) + bg(:,i) + gyro_noise(:,i); 
end

%% generate GPS data
lla0 = simData.lla0;
GPS_freq = 2; % Hz
ratio = ceil(IMU_freq / GPS_freq);
LLA = zeros(3, length(time));
for i = 1: length(time)
    LLA(:,i) = ned2lla(posI(:,i),lla0);
end
S_gps = diag([(1.35e-5/2);(1.63e-5/2);(2/2)]);
noisy_LLA = (LLA + S_gps * randn(size(LLA)));
GPS = nan(size(LLA));
for i = 1: ratio: length(time)
    GPS(:,i) = noisy_LLA(:,i); 
end
% GPS(:,1:3000) = nan(3,3000);
% GPS(:,6001:end) = nan(3,4000);
%% generate magnotometer data
% TODO need sensor datasheet
magI = [0.4891; 0.104; 0.866];
magB = zeros(size(bg));
sig_mag = 1e-8;
mag_noise = sig_mag * randn(size(magB));
for i = 1:size(magB,2)
    magB(:,i) = quaternion2DCM(q(:,i))' * magI ;
end
magB = magB + mag_noise;
%%
true_data.time = time;
true_data.acc = acc;
true_data.gyro = gyro;
true_data.ba = ba;
true_data.bg = bg;
true_data.Sa = Sa;
true_data.Sg = Sg;
true_data.LLA = LLA;
true_data.q = q;
true_data.GPS = GPS;
true_data.magB = magB;
save("sim_IMU_data.mat","true_data");
