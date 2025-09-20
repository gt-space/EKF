
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
N = length(t);
% DEFINE CONSTANTS %
we = 7.29211*10^(-5); % Earth's sidereal rotation rate [rad/s]
we = 0; % TODO only for current sim data
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

Q = dt*[nu_gv_mat, zeros(3,9);...
    zeros(3,3), nu_gu_mat, zeros(3,6);...
    zeros(3,6), nu_av_mat, zeros(3,3);...
    zeros(3,9), nu_au_mat];
Qq = dt*[nu_gv_mat zeros(3,3);
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
x_plus = x0;
P_plus = P0;
Pq_plus = Pq0;

Xplus_history = zeros(n,N);
Pplus_history = zeros(4,4,N);
% Run EKF
% running = true;
    
    % STEP 1: PROPAGATION
    % get sensor measurements at current time step
    %{
    a_meas = sens.sf_meas;
    w_meas = sens.omega_meas;
    lla_meas = sens.LLA_meas;
    %}
last_GPS_time = 0;
last_mag_time = 0;

for i = 1 : N
    % if i == sim_data.apogee_idx + 3000
    %     x_plus(11:13) = zeros(3,1);
    % end
    % iterates through each point
    a_meas = acc(:,i);
    w_meas = omega(:,i);
    w_hat = compute_what(x_plus,w_meas,we); % compute gyro estimate (accounting for sf & bias)
    ahat_n = compute_ahat(x_plus,a_meas); % compute accel estimate (accounting for sf, bias & frame transformation)
    Dt_GPS = dt*(i-1) - last_GPS_time;
    Dt_mag = dt*(i-1) - last_mag_time;
    [x_minus,P_minus,Pq_minus] = propagate(x_plus,P_plus,Pq_plus,w_hat,ahat_n,w_meas,a_meas,we,Q,Qq,dt);

    lla_meas = GPS(:,i);
    mag_meas = magB(:,i);
    press_meas = baro(i);

    if ~( sum(isnan(lla_meas)) || sum(isnan(mag_meas)) ) % both available
        [x_plus,P_plus,Pq_plus] = update_GPS_mag(x_minus,P_minus,Pq_minus,H,Hq,R,Rq,magI,lla_meas,mag_meas);
        last_GPS_time = dt * (i-1);
        last_mag_time = dt * (i-1);
    elseif ~sum(isnan(lla_meas))
        [x_plus,P_plus,Pq_plus] = update_GPS(x_minus,P_minus,Pq_minus,H,R,lla_meas);
        last_GPS_time = dt * (i-1);
    elseif ~sum(isnan(mag_meas))
        [x_plus,P_plus,Pq_plus] = update_mag(x_minus,P_minus,Pq_minus,Hq,Rq,magI,mag_meas);
        last_mag_time = dt * (i-1);
    else
        x_plus = x_minus;
        P_plus = P_minus;
        Pq_plus = Pq_minus;
    end
    % if i > 39060 && i < 44050
    %     x_plus(8:10) = x_minus(8:10); % not to mess up velocity
    % end
    % if Dt_GPS > 10 % lost GPS for more than 10 secs (arbitrary)
    %     x_plus(8:10) = x_minus(8:10); % not to mess up velocity
    % end
    if ~ sum(isnan(press_meas))
        [x_plus,P_plus,Pq_plus] = update_baro(x_plus,P_plus,Pq_plus,Hb,Rb,press_meas,hb);
    end
    if any(diag(P_plus)<0)
        fprintf("%d\n",i);
        P_plus = nearestPSD(P_plus); %+ diag(max(1e-6 - diag(P_plus), 0)); % heuristic to correct diagonal negative covariance
    end
    Xplus_history(:, i) = x_plus; % saves the value from xplus 1X22 vector for each dt
    Pplus_history(:,:, i) = P_plus(6:9,6:9);
end

att_error = zeros(1,N);
for k = 1:N
    q_true = sim_data.q(1:4,k);
    q_hat = Xplus_history(1:4,k);
    q_err = qmult(q_true,qconj(q_hat)); % right-invariant error quaternion
    att_error(k) = rad2deg(2*atan2(norm(q_err(2:4)),q_err(1)));  % error angle in radians
end
figure;
plot(t, att_error, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Attitude Error (degree)');
title('Attitude Error');
grid on;

ned_error = zeros(size(ned));
for i = 1:N
    ned_error(:,i) = lla2ned(Xplus_history(5:7,i),lla0)-ned(:,i);
end
figure
plot(t,ned_error(1:2,:), 'LineWidth', 2)
xlabel('time (s)')
ylabel('(m)')
title('Horizontal Position Error')
grid on;

figure
plot(t,-ned_error(3,:),'LineWidth', 2)
xlabel('time (s)')
ylabel('altitude error (m)')
title('Altitude Error')
grid on;

figure
plot(t,Xplus_history(8:10,:)-sim_data.velI(:,1:length(t)),'LineWidth', 2)
xlabel('time (s)')
ylabel('velocity error (m/s)')
title('Velocity Error')
legend('v1','v2','v3')
grid on;

figure
plot(t,squeeze(Pplus_history(1,1,:)),'LineWidth', 2)
hold on
plot(t,squeeze(Pplus_history(2,2,:)),'LineWidth', 2)
plot(t,squeeze(Pplus_history(3,3,:)),'LineWidth', 2)
plot(t,squeeze(Pplus_history(4,4,:)),'LineWidth', 2)
legend('P66','P77','P88','P99')
grid on;
hold off
%% Helper functions
function P_corrected = nearestPSD(P)
    [V, D] = eig((P + P')/2);           % Ensure symmetry and decompose
    D(D < 1e-10) = 1e-10;               % Clamp negative/near-zero eigenvalues
    P_corrected = V * D * V';              % Reconstruct matrix
end

function [x_plus,P_plus,Pq_plus] = update_baro(x_minus,P_minus,Pq_minus,Hb_func,Rb,press_meas,hb)
    % STEP 2: KALMAN GAIN
    Hb = Hb_func(x_minus);
    K = P_minus*Hb'/(Hb*P_minus*Hb' + Rb);

    % STEP 3: MEAS. UPDATE
    q_minus = x_minus(1:4);
    Delta_x = K*(press_meas - hb(x_minus)); % Step 3.1: compute DeltaX vector
    % q_plus = delta_q ⊗ q_minus
    q_plus = q_minus;
    p_plus = x_minus(5:7) + Delta_x(4:6);
    v_plus = x_minus(8:10) + Delta_x(7:9);
    bg_plus = x_minus(11:13);
    ba_plus = x_minus(14:16) + Delta_x(13:15);
    kg_plus = x_minus(17:19) + Delta_x(16:18);
    ka_plus = x_minus(20:22) + Delta_x(19:21);

    x_plus = [q_plus; p_plus; v_plus; bg_plus; ba_plus; kg_plus; ka_plus];
    P_plus = (eye(21,21) - K*Hb)*P_minus* (eye(21,21)-K * Hb)' + K*Rb*K';
    Pq_plus = Pq_minus;
end
function [x_plus,P_plus,Pq_plus] = update_mag(x_minus,P_minus,Pq_minus,Hq,Rq,magI,mag_meas)
    % STEP 2: KALMAN GAIN
    Kq = Pq_minus * Hq' / (Hq * Pq_minus * Hq' + Rq);

    % STEP 3: MEAS. UPDATE
    q_minus = x_minus(1:4);
    dv = magI - q2v(sandwich(q_minus,[0;mag_meas]));   
    err = Kq*dv;
    cq = err(1:3);
    cbg = err(4:6);
    dq = quaternion_exp(-cq/2);
    % q_plus = delta_q ⊗ q_minus
    q_plus = qmult(dq,q_minus);
    bg_plus = x_minus(11:13) - q2v(sandwich(qconj(q_plus),[0;cbg]));
    % if norm(bg_plus) > 1e-3
    %    bg_plus = zeros(3,1);
    % end
    x_plus = x_minus;
    x_plus(1:4) = q_plus;
    x_plus(11:13) = bg_plus;
    P_plus = P_minus;
    Pq_plus = (eye(6) - Kq * Hq) * Pq_minus * (eye(6) - Kq * Hq)' + Kq * Rq * Kq';   
end

function [x_plus,P_plus,Pq_plus] = update_GPS_mag(x_minus,P_minus,Pq_minus,H,Hq,R,Rq,magI,lla_meas,mag_meas)
    % STEP 2: KALMAN GAIN
    % underweighting
    if trace(P_minus(4:5,4:5)) + P_minus(6,6) > 1000
        beta = 0.25;
    else
        beta = 0;
    end
    W = (1+beta)*H*P_minus*H' + R;
    K = P_minus*H'/ W;
    % K = P_minus*H'/(H*P_minus*H' + R);
    Kq = Pq_minus * Hq' / (Hq * Pq_minus * Hq' + Rq);

    % STEP 3: MEAS. UPDATE
    q_minus = x_minus(1:4);
    Delta_x = K*(lla_meas - x_minus(5:7)); % Step 3.1: compute DeltaX vector
    dv = magI - q2v(sandwich(q_minus,[0;mag_meas]));   
    err = Kq*dv;
    cq = err(1:3);
    cbg = err(4:6);
    dq = quaternion_exp(-cq/2);
    % q_plus = delta_q ⊗ q_minus
    q_plus = qmult(dq,q_minus);
    p_plus = x_minus(5:7) + Delta_x(4:6);
    v_plus = x_minus(8:10) + Delta_x(7:9);
    % bg_plus = xk_minus(11:13) + Delta_x(10:12);
    bg_plus = x_minus(11:13) - q2v(sandwich(qconj(q_plus),[0;cbg]));
    ba_plus = x_minus(14:16) + Delta_x(13:15);
    kg_plus = x_minus(17:19) + Delta_x(16:18);
    ka_plus = x_minus(20:22) + Delta_x(19:21);
    % if norm(bg_plus) > 1e-3
    %    bg_plus = zeros(3,1);
    % end
    % if norm(kg_plus) > 0
    %    kg_plus = zeros(3,1);
    % end
    % if norm(ba_plus) > 0
    %    ba_plus = zeros(3,1);
    % end
    % if norm(ka_plus) > 0
    %    ka_plus = zeros(3,1);
    % end
    x_plus = [q_plus; p_plus; v_plus; bg_plus; ba_plus; kg_plus; ka_plus];
    P_plus = (eye(21,21) - K*H)*P_minus* (eye(21,21)-K * H)' + K*R*K';
    Pq_plus = (eye(6) - Kq * Hq) * Pq_minus * (eye(6) - Kq * Hq)' + Kq * Rq * Kq';   
end
function [x_minus,P_minus,Pq_minus] = propagate(x_plus,P_plus,Pq_plus,w_hat,ahat_n,w_meas,a_meas,we,Q,Qq,dt)
    % step 1.5: propagate to next step using Euler integration
    qdot = compute_qdot(x_plus,w_hat); % compute attitude derivative
    pdot = compute_lla_dot(x_plus); % compute position derivative
    vdot = compute_vdot(x_plus,ahat_n,we); % compute velocity derivative
    Pdot = compute_Pdot(x_plus,P_plus,Q,a_meas,w_meas,we); % compute covariance derivative
    Pqdot = compute_Pqdot(x_plus,Pq_plus,Qq,w_meas);
    [x_minus,P_minus,Pq_minus] = integrate(x_plus,P_plus,Pq_plus,dt,qdot,pdot,vdot,Pdot,Pqdot);
    % X0 = [x_plus;P_plus(:);Pq_plus(:)];
    % % RK4
    % [~,X] = ode45(@(t,x)dynamics(t,x,w_hat,ahat_n,w_meas,a_meas,we,Q,Qq),[0,dt],X0);
    % x_minus = X(end,1:22)';
    % P_minus = reshape(X(end,23:end-36)',21,21);
    % Pq_minus = reshape(X(end,end-35:end)',6,6);
end
function xdot = dynamics(t,x,w_hat,ahat_n,w_meas,a_meas,we,Q,Qq) % only used by RK4
    x_plus = x(1:22);
    P_plus = reshape(x(23:end-36),21,21);
    Pq_plus = reshape(x(end-35:end),6,6);
    qdot = compute_qdot(x_plus,w_hat); % compute attitude derivative
    pdot = compute_lla_dot(x_plus); % compute position derivative
    vdot = compute_vdot(x_plus,ahat_n,we); % compute velocity derivative
    Pdot = compute_Pdot(x_plus,P_plus,Q,a_meas,w_meas,we); % compute covariance derivative
    Pdot = Pdot(:);
    Pqdot = compute_Pqdot(x_plus,Pq_plus,Qq,w_meas);
    Pqdot = Pqdot(:);
    xdot = [qdot;pdot;vdot;zeros(12,1);Pdot;Pqdot];
end
function Pqdot = compute_Pqdot(x,Pq,Qq,w_meas)
    q = x(1:4);
    CB2I = quaternion2DCM(q);
    bg_minus = x(11:13);
    I_omega = CB2I*(w_meas - bg_minus);
    Fq = [zeros(3,3) -eye(3);
         zeros(3,3) skew(I_omega)];
    Gq = [CB2I zeros(3,3);
         zeros(3,3) -CB2I];
    Pqdot = Fq*Pq+ Pq*Fq' + Gq*Qq*Gq'; % compute covariance derivative
end

function wn = compute_wn(x,we)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); % north velocity
ve = x(9); % east velocity

[R_phi, R_lamb, ~, ~] = compute_radii(x);

vec1 = we*[cosd(phi); 0; -sind(phi)];
vec2 = [ve/(R_lamb + h); -vn/(R_phi + h); -(ve*tand(phi))/(R_lamb + h)];

wn = vec1 + vec2; % Eqn 7.70

end

function w_hat = compute_what(x,w_meas,we)

sf_g = x(17:19); % get gyro scale factor
bias_g = x(11:13); % get gyro bias
q = x(1:4); % get quaternion

D_bn = quaternion2DCM(q)'; % get DCM I2B from quaternion

w_n = compute_wn(x,we);

w_hat = 1./(1 + sf_g).*(w_meas - bias_g) - D_bn*w_n; % Eqn. 7.68b

end

function ahat_n = compute_ahat(x,a_meas)

sf_a = x(20:22); % get accel scale factor
bias_a = x(14:16); % get accel bias

ahat_b = 1./(1 + sf_a).*(a_meas - bias_a); % accel. estimate in body frame (eqn. 7.68j)

D_nb = quaternion2DCM(x(1:4)); % get DCM from body frame to NED frame 

ahat_n = D_nb*ahat_b; % get accel. estimate in NED frame
end

function qdot = compute_qdot(x,what)
q_bn = x(1:4); % get quaternion
w_quat = [0; what]; % convert angular rate into quaternion form
qdot = 1/2*qmult(q_bn,w_quat); % compute quaternion derivative
% qdot = 0.5 * what;
end


function lla_dot = compute_lla_dot(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); % north velocity
ve = x(9); % east velocity
vd = x(10); % down velocity

[R_phi, R_lamb, ~, ~] = compute_radii(x); % compute radii

phidot = vn/(R_phi + h); % Eqn 7.68c
lambdot = ve/((R_lamb + h)*cosd(phi)); % Eqn 7.68d
hdot = -vd; % Eqn 7.68e

lla_dot = [phidot; lambdot; hdot]; % assemble pdot vector
end

function vdot = compute_vdot(x,ahat_n,we)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); vd = x(10); % get vel components
an = ahat_n(1); ae = ahat_n(2); ad = ahat_n(3); % get accel components

[R_phi, R_lamb, ~, ~] = compute_radii(x); % compute radii
[ghat,~,~] = compute_g_dg(x); % get gravity estimate

vndot = -(ve/((R_lamb + h)*cosd(phi)) + 2*we)*ve*sind(phi) + (vn*vd)/(R_phi + h) + an; % Eqn 7.68f
vedot = (ve/((R_lamb + h)*cosd(phi)) + 2*we)*vn*sind(phi) + (ve*vd)/(R_lamb + h) + 2*we*vd*cosd(phi) + ae; % Eqn 7.68g
vddot = -ve^2/(R_lamb + h) - vn^2/(R_phi + h) - 2*we*ve*cosd(phi) + ghat + ad; % Eqn 7.68h
% vddot = -ve^2/(R_lamb + h) - vn^2/(R_phi + h) - 2*we*ve*cosd(phi) + ad; % Eqn 7.68h

vdot = [vndot; vedot; vddot]; % assemble vdot vector
end

function Pdot = compute_Pdot(x,P,Q,a_meas,w_meas,we)

F = compute_F(x,a_meas,w_meas,we);
G = compute_G(x);

Pdot = F*P + P*F' + G*Q*G';

end

function [g,dg_dphi,dg_dh] = compute_g_dg(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude
g = 9.780327*(1 + 5.3024*10^(-3)*sind(phi)^2 - 5.8*10^(-6)*sind(2*phi)^2) - ...
     (3.0877*10^(-6) - 4.4*10^(-9)*sind(phi)^2)*h + 7.2*10^(-14)*h^2; % Eqn 7.69c
% g = 9.81;

dg_dphi = 9.780327*(1.06048*10^(-2)*sind(phi)*cosd(phi)... % Eqn 7.84a
    - 4.64*10^(-5)*(sind(phi)*cosd(phi)^3 - sind(phi)^3*cosd(phi)))...
    + 8.8*10^(-9)*h*sind(phi)*cosd(phi);

dg_dh = -3.0877*10^(-6) + 4.4*10^(-9)*sind(phi)^2 + 1.44*10^(-13)*h; % Eqn 7.84b

end

function [R_phi, R_lamb, dR_phi_dphi, dR_lamb_dphi] = compute_radii(x)

% WGS84 model parameters (source: https://en.wikipedia.org/wiki/World_Geodetic_System)
a = 6378137; % semi-major axis [m]
b = 6356752.31425; % semi-minor axis [m]
e2 = 1 - (b/a)^2; % eccentricity

phi = x(5); % get latitude

num1 = a*(1-e2);
den1 = (1-e2*sind(phi)^2)^(3/2);

den2 = (1-e2*sind(phi)^2)^(1/2);

R_phi = num1/den1; % Eqn 7.69a
R_lamb = a/den2; % Eqn 7.69b

num3 = 3*a*(1-e2)*e2*sind(phi)*cosd(phi);
den3 = (1-e2*sind(phi)^2)^(5/2);

num4 = a*e2*sind(phi)*cosd(phi);
den4 = (1-e2*sind(phi)^2)^(3/2);

dR_lamb_dphi = num4/den4; % Eqn 7.75a
dR_phi_dphi = num3/den3; % Eqn 7.75b

end

function F = compute_F(x,a_meas,w_meas,we)

sf_g = x(17:19); % get gyro scale factor
sf_a = x(20:22); % get accel scale factor
bias_g = x(11:13); % get gyro bias
bias_a = x(14:16); % get accel bias
q = x(1:4); % get quaternion

D_nb = quaternion2DCM(q); % get DCM from body frame to NED frame
D_bn = D_nb'; 

F11_vec = -1./(1 + sf_g).*(w_meas - bias_g);
F11 = skew(F11_vec); % Eqn. 7.79a

dwdp = compute_dwdp(x,we);
F12 = -D_bn*dwdp; % Eqn. 7.79a

dwdv = compute_dwdv(x);
F13 = -D_bn*dwdv; % Eqn. 7.79b

F14 = diag(-1./(1 + sf_g)); % Eqn. 7.79b

Omega = diag(w_meas);
Bg = diag(bias_g);
F16 = -(Omega - Bg); % Eqn. 7.79b

F22 = compute_dpdot_dp(x); % Eqn. 7.79c
F23 = compute_dpdot_dv(x); % Eqn. 7.79c

ahat_n = compute_ahat(x,a_meas);
ahat_b = D_bn*ahat_n;
F31 = -D_nb*skew(ahat_b); % Eqn. 7.79d

F32 = compute_dvdot_dp(x,we); % Eqn. 7.79d
F33 = compute_dvdot_dv(x,we); % Eqn. 7.79d

F35 = -D_nb*diag(1./(1 + sf_a)); % Eqn. 7.79e
F37 = -D_nb*(diag(a_meas) - diag(bias_a)); % Eqn. 7.79e

F = [F11, F12, F13, F14, zeros(3,3), F16, zeros(3,3);... % Eqn 7.78a
    zeros(3,3), F22, F23, zeros(3,12);...
    F31, F32, F33, zeros(3,3), F35, zeros(3,3), F37;...
    zeros(12,21)];
end

function G = compute_G(x)

sf_g = x(17:19); % get gyro scale factor
sf_a = x(20:22); % get accel scale factor
q = x(1:4); % get quaternion

D_bn = quaternion2DCM(q)'; % get DCM from quaternion
D_nb = D_bn'; % get DCM from body frame to NED frame

G11 = -diag(1./(1 + sf_g));
G33 = -D_nb*diag(1./(1 + sf_a));

G = [G11, zeros(3,9);... % Eqn 7.78b
    zeros(3,12);...
    zeros(3,6), G33, zeros(3,3);...
    zeros(3,3), eye(3), zeros(3,6);...
    zeros(3,9), eye(3);...
    zeros(6,12)];

end

function dwdp = compute_dwdp(x,we)
%Eqn 7.74a

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); % get vel components

[R_phi, R_lamb, dR_phi_dphi, dR_lamb_dphi] = compute_radii(x);

% Eqn 7.74a
m11 = -we*sind(phi) - ve/(R_lamb + h)^2 * dR_lamb_dphi;
m13 = -ve/(R_lamb + h)^2;
m21 = vn/(R_phi + h)^2 * dR_phi_dphi;
m23 = vn/(R_phi + h)^2;
m31 = -we*cosd(phi)-(ve*secd(phi)^2)/(R_lamb + h) + (ve*tand(phi))/(R_lamb + h)^2*dR_lamb_dphi;
m33 = (ve*tand(phi))/(R_lamb + h)^2;

% assemble dwdp matrix
dwdp = [m11, 0, m13;...
    m21, 0, m23;...
    m31, 0, m33];

end

function dwdv = compute_dwdv(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude

[R_phi, R_lamb, ~, ~] = compute_radii(x);

% Eqn 7.74b
m12 = 1/(R_lamb + h);
m21 = -1/(R_phi + h);
m32 = -tand(phi)/(R_lamb + h);

% assemble dwdv matrix
dwdv = [0, m12, 0;...
    m21, 0, 0;...
    0, m32, 0];

end

function dpdot_dp = compute_dpdot_dp(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); % get vel components

[R_phi, R_lamb, dR_phi_dphi, dR_lamb_dphi] = compute_radii(x);

% Eqn 7.80a
m11 = -vn/(R_phi + h)^2*dR_phi_dphi;
m13 = -vn/(R_phi + h)^2;
m21 = -(ve*secd(phi))/(R_lamb + h)^2*dR_lamb_dphi + (ve*secd(phi)*tand(phi))/(R_lamb + h);
m23 = -ve*secd(phi)/(R_lamb + h)^2;

% assemble dpdot_dp matrix
dpdot_dp = [m11, 0, m13;...
    m21, 0, m23;...
    0, 0, 0];

end

function dpdot_dv = compute_dpdot_dv(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude

[R_phi, R_lamb, ~, ~] = compute_radii(x);

% Eqn 7.80b
m11 = 1/(R_phi + h);
m22 = secd(phi)/(R_lamb + h);
m33 = -1;

% assemble dpdot_dv matrix
dpdot_dv = diag([m11, m22, m33]);
end

function dvdot_dp = compute_dvdot_dp(x,we)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); vd = x(10); % get vel components

[R_phi, R_lamb, dR_phi_dphi, dR_lamb_dphi] = compute_radii(x);

[~,dg_dphi,dg_dh] = compute_g_dg(x);

Y11 = -(ve^2*secd(phi)^2)/(R_lamb + h) + (ve^2*tand(phi))/(R_lamb + h)^2*dR_lamb_dphi ... % Eqn 7.82a
    - 2*we*ve*cosd(phi) - (vn*vd)/(R_phi + h)^2 * dR_phi_dphi;
Y13 = (ve^2*tand(phi))/(R_lamb + h)^2 - (vn*vd)/(R_phi + h)^2; % Eqn 7.82b
Y21 = (ve*vn*secd(phi)^2)/(R_lamb + h) - (ve*vn*tand(phi))/(R_lamb + h)^2*dR_lamb_dphi + ... % Eqn 7.82c
    2*we*vn*cosd(phi) - (ve*vd)/(R_lamb + h)^2*dR_lamb_dphi - 2*we*vd*sind(phi);
Y23 = -ve*((vn*tand(phi) + vd)/(R_lamb + h)^2); % Eqn 7.82d
Y31 = ve^2/(R_lamb + h)^2*dR_lamb_dphi + vn^2/(R_phi + h)^2*dR_phi_dphi + 2*we*ve*sind(phi) + dg_dphi; % Eqn 7.82e
Y33 = ve^2/(R_lamb + h)^2 + vn^2/(R_phi + h)^2 + dg_dh; % Eqn 7.82f

% assemble dvdot_dp matrix
dvdot_dp = [Y11, 0, Y13;...
    Y21, 0, Y23;...
    Y31, 0, Y33];
end

function dvdot_dv = compute_dvdot_dv(x,we)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); vd = x(10); % get vel components

[R_phi, R_lamb,~,~] = compute_radii(x);

Z11 = vd/(R_phi + h); % Eqn 7.83a
Z12 = (-2*ve*tand(phi))/(R_lamb + h) + 2*we*sind(phi); % Eqn 7.83a
Z13 = vn/(R_phi + h); % Eqn 7.83a
Z21 = (ve*tand(phi))/(R_lamb + h) + 2*we*sind(phi); % Eqn 7.83b
Z22 = (vd + vn*tand(phi))/(R_lamb + h); % Eqn 7.83b
Z23 = ve/(R_lamb + h) + 2*we*cosd(phi); % Eqn 7.83b
Z31 = (-2*vn)/(R_phi + h); % Eqn 7.83c
Z32 = (-2*ve)/(R_lamb + h) - 2*we*cosd(phi); % Eqn 7.83c

% assemble dvdot_dv matrix
dvdot_dv = [Z11, Z12, Z13;...
    Z21, Z22, Z23;...
    Z31, Z32,  0];

end


function [x_minus,P_minus,Pq_minus] = integrate(x,P,Pq,dt,qdot,pdot,vdot,Pdot,Pqdot)

% assemble xdot vector
xdot = [qdot; pdot; vdot; zeros(12,1)];

x_minus = x + dt*xdot; % use Euler integration
q_next = x_minus(1:4); % get propagated quaternion

q_norm = norm(q_next);
if q_norm > 1.02 % if quaternion norm is above 1% greater than unity
%     warning('Quaternion norm of %.2f is > 2 percent of unity, consider decreasing integrator step size.',norm(q_next))
end
x_minus(1:4) = q_next/q_norm; % normalize quaternion

P_minus = P + dt*Pdot; % use Euler integration
Pq_minus = Pq + dt*Pqdot; % use Euler integration

% xdot = [pdot; vdot; zeros(12,1)];
% q = x(1:4);
% q_next = qmult(q, quaternion_exp(qdot)); % use exponential map to get unit quaternion
% q_norm = norm(q_next);
% if q_norm > 1.01 % if quaternion norm is above 1% greater than unity
% %     warning('Quaternion norm of %.2f is > 1 percent of unity, consider decreasing integrator step size.',norm(q_next))
% end
% q_minus = q_next/q_norm; % normalize quaternion
% 
% x_minus = [q_minus; x(5:end) + dt*xdot];
% P_minus = P + dt*Pdot; % use Euler integration
% Pq_minus = Pq + dt*Pqdot; % use Euler integration

end

function Xi = createXiMat(q)
% create Xi matrix from scalar-first quaternion
% used eqn 2.88 from sat ADCS textbook and adapted it for scalar-first quaternion

q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);

Xi = [q1, -q4, q3; q4, q1, -q2; -q3, q2, q1; -q2, -q3, -q4];

end