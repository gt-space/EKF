
% CODE FOR VALIDATING ALGORITHIM WITH SAME SAMPLE RATE FOR IMU & GPS
close all; clear; clc;

load('sim_data_4_EKF_testing.mat');
acc = sim_data.acc;
omega = sim_data.gyro;
lla = sim_data.LLA;
GPS = sim_data.GPS;
magB = sim_data.magB;
magI = sim_data.magI;
t = sim_data.time;
N = length(t);
% DEFINE CONSTANTS %
we = 7.29211*10^(-5); % Earth's sidereal rotation rate [rad/s]

% define H
H = zeros(3,21); H(:,4:6) = eye(3); % Eqn for Hk in Table 7.2

Hq = [skew(magI) zeros(3,3)];
% define R
R = diag([1e-4 1e-4 2]); % MODIFY TO MODEL GNSS NOISE
Rq = diag([1e-7 1e-7 1e-7]);
% define Q ( eqn 7.76b)
% nu_gv_mat = 6.8256*10^(-10)*eye(3); % using rad/s
% nu_gu_mat = 3.1002*10^(-7)*eye(3); % using rad/s
% nu_av_mat = 1.6867*10^(-4)*eye(3); % using m/s^2
% nu_au_mat = 1.3058*10^(-6)*eye(3); % using m/s^2

% MODIFY TO MODEL GYRO/ACCEL NOISE
nu_gv_mat = 0.1*eye(3); % using rad/s
nu_gu_mat = 10^(-6)*eye(3); % using rad/s
nu_av_mat = 0.1*eye(3); % using m/s^2
nu_au_mat = 10^(-6)*eye(3); % using m/s^2

Q = [nu_gv_mat, zeros(3,9);...
    zeros(3,3), nu_gu_mat, zeros(3,6);...
    zeros(3,6), nu_av_mat, zeros(3,3);...
    zeros(3,9), nu_au_mat];
Qq = [nu_gv_mat zeros(3,3);
      zeros(3,3) nu_gu_mat];
Iq = eye(6);
% DEFINE PARAMETERS
sr = 100; % save sample rate variable
dt = 1/sr; % get time step

% DEFINE INITIAL CONDITIONS %
% lla0 = [30.9275; -81.51472222222; 45];
lla0 = lla(:,1);
q0 = [0.707106781186548; 0; 0.707106781186547; 0];

% attitude, pos, vel, g bias, a bias, g sf, a sf
x0 = [q0; lla0; zeros(3,1); zeros(3,1); zeros(3,1); zeros(3,1); zeros(3,1)];
n = size(x0,1);

att_unc0 = 0.01; % initial attitude uncertainty [rad]
pos_unc0 = 0.01; % initial position uncertainty [m]
vel_unc0 = 0.001; % initial velocity uncertainty [m/s]
gbias_unc0 = 0.0001; % initial gyro bias uncertainty [rad/s]
abias_unc0 = 0.0001; % initial accel bias uncertainty [rad/s]
gsf_unc0 = 0.0001; % gyro scale factor uncertainty [rad/s]
asf_unc0 = 0.005; % accel scale factor uncertainty [rad/s]

% create initial covariance matrix
P0 = diag([att_unc0*ones(1,3), pos_unc0*ones(1,3), vel_unc0*ones(1,3), gbias_unc0*ones(1,3),...
    abias_unc0*ones(1,3),gsf_unc0*ones(1,3), asf_unc0*ones(1,3)]);
Pq0 = diag([att_unc0*ones(1,3),abias_unc0*ones(1,3)]);
% set up initial conditions for the filter
Pk_plus = P0;
Pq_plus = Pq0;
xk_plus = x0;

Xplus_history = zeros(n,N);
% Run EKF
% running = true;
    
    % STEP 1: PROPAGATION
    % get sensor measurements at current time step
    %{
    a_meas = sens.sf_meas;
    w_meas = sens.omega_meas;
    lla_meas = sens.LLA_meas;
    %}
for i = 1 : N
    % iterates through each point
    a_meas = acc(:,i);
    w_meas = omega(:,i);
    lla_meas = lla(:,i);
    
    w_hat = compute_what(xk_plus,w_meas,we); % compute gyro estimate (accounting for sf & bias)
    qdot = compute_qdot(xk_plus,w_hat); % compute attitude derivative
    pdot = compute_pdot(xk_plus); % compute position derivative
    
    ahat_n = compute_ahat(xk_plus,a_meas); % compute accel estimate (accounting for sf, bias & frame transformation)
    vdot = compute_vdot(xk_plus,ahat_n,we); % compute velocity derivative
    Pdot = compute_Pdot(xk_plus,Pk_plus,Q,a_meas,w_meas,we); % compute covariance derivative
    
    % step 1.5: propagate to next step using Euler integration
    [xk_minus,Pk_minus,q_norm] = integrate(xk_plus,Pk_plus,dt,qdot,pdot,vdot,Pdot);

    q = xk_plus(1:4);
    CB2I = quaternion2DCM(q);
    bg_minus = xk_minus(11:13);
    I_omega = CB2I*(w_meas - bg_minus);
    Fq = [zeros(3,3) -eye(3);
         zeros(3,3) skew(I_omega)];
    Gq = [CB2I zeros(3,3);
         zeros(3,3) -CB2I];
    Pqdot = Fq*Pq_plus+ Pq_plus*Fq' + Gq*Qq*Gq'; % compute covariance derivative
    Pq_minus = Pq_plus + dt * Pqdot;
    % STEP 2: KALMAN GAIN
    Kk = Pk_minus*H'/(H*Pk_minus*H' + R);
    Kq = Pq_minus * Hq' / (Hq * Pq_minus * Hq' + Rq);

    % STEP 3: MEAS. UPDATE
    q_minus = xk_minus(1:4);
    Delta_x = Kk*(lla_meas - xk_minus(5:7)); % Step 3.1: compute DeltaX vector
    dv = magI - q2v(sandwich(q_minus,[0;magB(:,i)]));   
    err = Kq*dv;
    cq = err(1:3);
    cbg = err(4:6);
    dq = quaternion_exp(-cq/2);
    % q_plus = delta_q ⊗ q_minus
    q_plus = qmult(dq,q_minus);
    p_plus = xk_minus(5:7) + Delta_x(4:6);
    v_plus = xk_minus(8:10) + Delta_x(7:9);
    % bg_plus = xk_minus(11:13) + Delta_x(10:12);
    bg_plus = xk_minus(11:13) - q2v(sandwich(qconj(q_plus),[0;cbg]));
    ba_plus = xk_minus(14:16) + Delta_x(13:15);
    kg_plus = xk_minus(17:19) + Delta_x(16:18);
    ka_plus = xk_minus(20:22) + Delta_x(19:21);
    
    xk_plus = [q_plus; p_plus; v_plus; bg_plus; ba_plus; kg_plus; ka_plus];
    Xplus_history(:, i) = xk_plus; % saves the value from xplus 1X22 vector for each dt
    Pk_plus = (eye(21,21) - Kk*H)*Pk_minus;
    Pq_plus = (Iq - Kq * Hq) * Pq_minus * (Iq - Kq * Hq)' + Kq * Rq * Kq';
end
att_error = zeros(1,N);
for k = 1:N
    q_true = sim_data.q(1:4,k);
    q_hat = Xplus_history(1:4,k);
    q_err = qmult(q_true,qconj(q_hat)); % right-invariant error quaternion
    att_error(k) = 2*acos(abs(q_err(1)));  % error angle in radians
end
figure;
plot(t, att_error, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Attitude Error (rad)');
title('Attitude Error');
grid on;

figure
plot(t,Xplus_history(5:6,:)-lla(1:2,:), 'LineWidth', 2)
xlabel('time (s)')
ylabel('degree')
title('Horizontal Position Error')
grid on;

figure
plot(t,Xplus_history(7,:)-lla(3,:),'LineWidth', 2)
xlabel('time (s)')
ylabel('altitude error (m)')
title('Altitude Error')
grid on;

%% Helper functions
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

D_bn = quaternion2DCM(q)'; % get DCM from quaternion
sf_mat_g = diag(sf_g); % create scale factor matrix for gyro

w_n = compute_wn(x,we);
% w_n = zeros(3,1); 

w_hat = (eye(3) - sf_mat_g)*(w_meas - bias_g) - D_bn*w_n; % Eqn. 7.68b

end

function ahat_n = compute_ahat(x,a_meas)

sf_a = x(20:22); % get accel scale factor
bias_a = x(14:16); % get accel bias

sf_mat_a = diag(sf_a); % create scale factor matrix for gyro

ahat_b = (eye(3) - sf_mat_a)*(a_meas - bias_a); % accel. estimate in body frame (eqn. 7.68j)

D_bn = quaternion2DCM(x(1:4))'; % get DCM from quaternion
D_nb = D_bn'; % get DCM from body frame to NED frame

ahat_n = D_nb*ahat_b; % get accel. estimate in NED frame
end

function qdot = compute_qdot(x,what)
q_bn = x(1:4); % get quaternion
w_quat = [0; what]; % convert angular rate into quaternion form
qdot = 1/2*qmult(q_bn,w_quat); % compute quaternion derivative
end


function pdot = compute_pdot(x)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); % north velocity
ve = x(9); % east velocity
vd = x(10); % down velocity

[R_phi, R_lamb, ~, ~] = compute_radii(x); % compute radii

phidot = vn/(R_phi + h); % Eqn 7.68c
lambdot = ve/((R_lamb + h)*cosd(phi)); % Eqn 7.68d
hdot = -vd; % Eqn 7.68e

pdot = [phidot; lambdot; hdot]; % assemble pdot vector
end

function vdot = compute_vdot(x,ahat_n,we)

phi = x(5); % get latitude
h = x(7); % get height/altitude
vn = x(8); ve = x(9); vd = x(10); % get vel components
an = ahat_n(1); ae = ahat_n(2); ad = ahat_n(3); % get accel components

[R_phi, R_lamb, ~, ~] = compute_radii(x); % compute radii
[ghat,~,~] = compute_g(x); % get gravity estimate

vndot = -((ve/(R_lamb + h)*cosd(phi)) + 2*we)*ve*sind(phi) + (vn*vd)/(R_phi + h) + an; % Eqn 7.68f
vedot = ((ve/(R_lamb + h)*cosd(phi)) + 2*we)*vn*sind(phi) + (ve*vd)/(R_lamb + h) + 2*we*vd*cosd(phi) + ae; % Eqn 7.68g
vddot = -ve^2/(R_lamb + h) - vn^2/(R_phi + h) - 2*we*ve*cosd(phi) + ghat + ad; % Eqn 7.68h

vdot = [vndot; vedot; vddot]; % assemble vdot vector
end

function Pdot = compute_Pdot(x,P,Q,a_meas,w_meas,we)

F = compute_F(x,a_meas,w_meas,we);
G = compute_G(x);

Pdot = F*P + P*F' + G*Q*G';

end

function [g,dg_dphi,dg_dh] = compute_g(x)

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
e = sqrt(1 - b^2/a^2); % eccentricity

phi = x(5); % get latitude

num1 = a*(1-e^2);
den1 = (1-e^2*sind(phi)^2)^(3/2);

den2 = (1-e^2*sind(phi)^2)^(1/2);

R_phi = num1/den1; % Eqn 7.69a
R_lamb = a/den2; % Eqn 7.69b

num3 = 3*a*(1-e^2)*e^2*sind(phi)*cosd(phi);
den3 = (1-e^2*sind(phi)^2)^(5/2);

num4 = a*e^2*sind(phi)*cosd(phi);
den4 = (1-e^2*sind(phi)^2)^(3/2);

dR_lamb_dphi = num4/den4; % Eqn 7.75a
dR_phi_dphi = num3/den3; % Eqn 7.75b

end

function F = compute_F(x,a_meas,w_meas,we)

sf_g = x(17:19); % get gyro scale factor
sf_a = x(20:22); % get accel scale factor
bias_g = x(11:13); % get gyro bias
bias_a = x(14:16); % get accel bias
q = x(1:4); % get quaternion

D_bn = quaternion2DCM(q)'; % get DCM from quaternion
D_nb = D_bn'; % get DCM from body frame to NED frame
sf_mat_g = diag(sf_g); % create scale factor matrix for gyro
sf_mat_a = diag(sf_a); % create scale factor matrix for accel

F11_vec = -(eye(3) - sf_mat_g)*(w_meas - bias_g);
F11 = createSkewSymmMat(F11_vec); % Eqn. 7.79a

dwdp = compute_dwdp(x,we);
F12 = -D_bn*dwdp; % Eqn. 7.79a

dwdv = compute_dwdv(x);
F13 = -D_bn*dwdv; % Eqn. 7.79b

F14 = -(eye(3) - sf_mat_g); % Eqn. 7.79b

Omega = diag(w_meas);
Bg = diag(bias_g);
F16 = -(Omega - Bg); % Eqn. 7.79b

F22 = compute_dpdot_dp(x); % Eqn. 7.79c
F23 = compute_dpdot_dv(x); % Eqn. 7.79c

ahat_n = compute_ahat(x,a_meas);
ahat_b = D_bn*ahat_n;
F31 = -D_nb*createSkewSymmMat(ahat_b); % Eqn. 7.79d

F32 = compute_dvdot_dp(x,we); % Eqn. 7.79d
F33 = compute_dvdot_dv(x,we); % Eqn. 7.79d

F35 = -D_nb*(eye(3) - sf_mat_a); % Eqn. 7.79e
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
sf_mat_g = diag(sf_g); % create scale factor matrix for gyro
sf_mat_a = diag(sf_a); % create scale factor matrix for accel

G11 = -(eye(3) - sf_mat_g);
G33 = -D_nb*(eye(3) - sf_mat_a);

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

[~,dg_dphi,dg_dh] = compute_g(x);

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

function ssMat = createSkewSymmMat(vec)
%input is a 3x1 vector
%output is a skew symmetric matrix using that vector

ssMat = [0, -vec(3), vec(2);
    vec(3), 0, -vec(1);
    -vec(2), vec(1), 0];
end

function [x_next,P_next,q_norm] = integrate(x,P,dt,qdot,pdot,vdot,Pdot)

% assemble xdot vector
xdot = [qdot; pdot; vdot; zeros(12,1)];

x_next = x + dt*xdot; % use Euler integration
q_next = x_next(1:4); % get propagated quaternion

q_norm = norm(q_next);
if q_norm > 1.02 % if quaternion norm is above 1% greater than unity
%     warning('Quaternion norm of %.2f is > 2 percent of unity, consider decreasing integrator step size.',norm(q_next))
end

x_next(1:4) = q_next/q_norm; % normalize quaternion

P_next = P + dt*Pdot; % use Euler integration
end

function Xi = createXiMat(q)
% create Xi matrix from scalar-first quaternion
% used eqn 2.88 from sat ADCS textbook and adapted it for scalar-first quaternion

q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);

Xi = [q1, -q4, q3; q4, q1, -q2; -q3, q2, q1; -q2, -q3, -q4];

end