
% CODE FOR VALIDATING ALGORITHIM WITH SAME SAMPLE RATE FOR IMU & GPS
close all; clear;

load('sim_data_4_EKF_testing.mat');
acc = sim_data.acc;
omega = sim_data.gyro;
ned = sim_data.NED;
GPS = sim_data.GPS;

magB = sim_data.magB;
magI = sim_data.magI;
baro = sim_data.baro;
t = sim_data.time;
N = floor(length(t)/1);

% DEFINE CONSTANTS %
we = 7.29211*10^(-5); % Earth's sidereal rotation rate [rad/s]
% DEFINE PARAMETERS
IMU_freq = 100; % save sample rate variable
dt = 1/IMU_freq; % get time step
GPS_freq = 1;
mag_freq = 10;
baro_freq = 20;

% define H
H = zeros(3,21); H(:,4:6) = eye(3); % Eqn for Hk in Table 7.2
p0 = 102715.47296217596;
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

hb = @(x) pressure_function(x(7),p0,a);
Hb = @(x) [zeros(1,5) pressure_derivative(x(7),p0,a) zeros(1,15)]; % barometer measurement function Jacobian
Hq = [skew(magI) zeros(3,3)]; % magnetometer measurement Jacobian and gyro bias
% define R

R = diag([1.35e-5;1.65e-5;2]); % MODIFY TO MODEL GNSS NOISE
Rq = diag([3.2e-7 4.1e-7 3.2e-7]);
Rb = 2.5e-3; % barometer noise variance


% define Q ( eqn 7.76b)
% nu_gv_mat = 6.8256*10^(-10)*eye(3); % using rad/s
% nu_gu_mat = 3.1002*10^(-7)*eye(3); % using rad/s
% nu_av_mat = 1.6867*10^(-4)*eye(3); % using m/s^2
% nu_au_mat = 1.3058*10^(-6)*eye(3); % using m/s^2
% MODIFY TO MODEL GYRO/ACCEL NOISE
nu_gv_mat = deg2rad(12e-3)^1*eye(3); % using rad/s
nu_gu_mat = deg2rad(3/3600)^1*eye(3); % using rad/s
nu_av_mat = (200e-6*9.81)^1*eye(3); % using m/s^2
nu_au_mat = (40e-6 *9.8)^1*eye(3); % using m/s^2
Q = 10*dt*[nu_gv_mat, zeros(3,9);...
    zeros(3,3), nu_gu_mat, zeros(3,6);...
    zeros(3,6), nu_av_mat, zeros(3,3);...
    zeros(3,9), nu_au_mat];
Qq = 10*dt*[nu_gv_mat zeros(3,3);
      zeros(3,3) nu_gu_mat];

% DEFINE INITIAL CONDITIONS %
% lla0 = [30.9275; -81.51472222222; 45];
lla0 = sim_data.LLA0;
q0 = [0.707106781186548; 0; 0.707106781186547; 0];
% attitude, pos, vel, g bias, a bias, g sf, a sf
x0 = [q0; lla0; zeros(3,1); zeros(3,1); zeros(3,1); zeros(3,1); zeros(3,1)];
n = size(x0,1);

att_unc0 = 1e-4; % initial attitude uncertainty [rad]
pos_unc0 = 1e-4; % initial position uncertainty [m]
vel_unc0 = 1e-4; % initial velocity uncertainty [m/s]
gbias_unc0 = 1e-4; % initial gyro bias uncertainty [rad/s]
abias_unc0 = 1e-4; % initial accel bias uncertainty [rad/s]
gsf_unc0 = 1e-4; % gyro scale factor uncertainty [rad/s]
asf_unc0 = 1e-4; % accel scale factor uncertainty [rad/s]

% create initial covariance matrix
P0 = diag([att_unc0*ones(1,3), pos_unc0*ones(1,3), vel_unc0*ones(1,3), gbias_unc0*ones(1,3),...
    abias_unc0*ones(1,3),gsf_unc0*ones(1,3), asf_unc0*ones(1,3)]);
Pq0 = diag([att_unc0*ones(1,3),gbias_unc0*ones(1,3)]);
% set up initial conditions for the filter

x_plus_history = single(zeros(n,N));
x_plus_history(:,1) = single(x0);
P_prev = single(P0);
Pq_prev = single(Pq0);
drogue_flag = 1;

%% EKF
for i = 1:N
    % TODO add Euler rate to angular velocity conversion
    a_meas = acc(:,i);
    w_meas = omega(:,i);
    lla_meas = GPS(:,i);
    mag_meas = magB(:,i);
    press_meas = baro(i);
    [x_plus,P_plus,Pq_plus,deploy_drogue,deploy_chute] = ...
        EKF_update(x_plus_history(:,i),P_prev,Pq_prev,Q,Qq,H,Hq,hb,Hb,R,Rq,Rb,a_meas,w_meas,lla_meas,mag_meas,press_meas,we,magI,dt);
    x_plus_history(:,i+1) = x_plus;

    if deploy_drogue == 1 && drogue_flag == 1
        fprintf("%d\n",i);
        fprintf("Deplying drogue at altitude %d, down speed %d.\n",x_plus(7),x_plus(10));
        drogue_flag = 0;
    end
    if mod(i,1000000) == 0 % periodically check if P is PSD
        P_plus = nearestPSD(P_plus);
    end

    P_prev = P_plus;
    Pq_prev = Pq_plus;
end

%% Plotting
att_error = zeros(1,N);
for k = 1:N
    q_true = sim_data.q(1:4,k);
    q_hat = x_plus_history(1:4,k);
    q_err = qmult(q_true,qconj(q_hat)); % right-invariant error quaternion
    att_error(k) = rad2deg(2*atan2(norm(q_err(2:4)),q_err(1)));  % error angle in radians
end

figure;
plot(t(1:N), att_error, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Attitude Error (degree)');
title('Attitude Error');
grid on;

ned_error = zeros(size(ned,1),N);
for i = 1:N
    ned_error(:,i) = lla2ned(x_plus_history(5:7,i),lla0)-ned(:,i);
end
figure
plot(t(1:N),ned_error(1:2,:), 'LineWidth', 2)
xlabel('time (s)')
ylabel('(m)')
title('Horizontal Position Error')
grid on;

figure
plot(t(1:N),-ned_error(3,:),'LineWidth', 2)
xlabel('time (s)')
ylabel('altitude error (m)')
title('Altitude Error')
grid on;

figure
plot(t(1:N),x_plus_history(8:10,2:end)-sim_data.velI(:,1:N),'LineWidth', 2)
xlabel('time (s)')
ylabel('velocity error (m/s)')
title('Velocity Error')
legend('v1','v2','v3')
grid on;

figure
plot(t(1:N),x_plus_history(10,2:end)-sim_data.velI(3,1:N),'LineWidth', 2)
xlabel('time (s)')
ylabel('velocity error (m/s)')
title('Vertical Velocity Error')
grid on;
