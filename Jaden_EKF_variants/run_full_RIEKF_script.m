clear;
load("sim_data_4_EKF_testing.mat")
measurements = [true_data.magB;true_data.GPS];
% IMU_freq = 1000; % Hz
IMU_freq = 100; % Hz
GPS_freq = 25; % Hz
Baro_freq = 100; % Hz
Mag_freq = 100; % Hz
dt = 1/IMU_freq;
ratio = IMU_freq / GPS_freq;

time = double(true_data.time);
N = length(time);
magI = [0.4891; 0.104; 0.866];
get_inert_mag_vec = @(x) magI; % assume constant mag vector
U = [true_data.gyro;true_data.acc];

hq = @(x) quaternion2DCM(x(1:4))'*magI ; % magnetometer measures magB
h = @(rs) rs(1:3); % rs is the remaining state, measures LLA position
hb = @(rs) rs(3); % barometer TODO
f = @(t,x,uB,IMU) qrv_dynamics(t,x,uB,IMU);
Dh = @(rs) h_jacobian(rs);
Df = @(x,u) f_jacobian(x,u);
q0 = true_data.q(:,1);
r0 = true_data.LLA(:,1);
v0 = [0;0;0];
x0 = [q0;r0;v0];

n = size(x0,1);
nr = n - 4;
Pq0 = blkdiag(1e-8*eye(3), 1e-6*eye(3)); % min error state dim
Pr0 = blkdiag(1e-4*eye(3), 1e-6*eye(3), 1e-6*eye(3),1e-6*eye(3));
bg0 = [0;0;0];
ba0 = [0;0;0];
Sa0 = [1;1;1];
IMU0 = [bg0;ba0;Sa0];

%%
get_Q = @(x) compute_Q(x,dt);
xhat0 = x0 + 1e-5*randn(size(x0));
xhat0(1:4) = xhat0(1:4)/norm(xhat0(1:4));
% xhat0 = x0;
[Xhat,IMU,Pq,Pr] = hybrid_full_RIEKF_offline(f,Df,h,Dh,xhat0,IMU0,Pq0,Pr0,time,U,...
                                measurements,get_Q,@get_R,get_inert_mag_vec,ratio);
att_error = zeros(length(time),1);
position_error = Xhat(5:7,:) - true_data.LLA;
for k = 1:length(time)
    q_true = true_data.q(1:4,k);
    q_hat = Xhat(1:4,k);
    q_err = qmult(qconj(q_hat), q_true); % right-invariant error quaternion
    q_err = q_err / norm(q_err);
    att_error(k) = 2*acos(min(1,abs(q_err(1))));  % error angle in radians
end
% Plot attitude (quaternion) error
figure;
plot(time, att_error, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Attitude Error (rad)');
title('RIEKF: Attitude Error');
grid on;
figure;
plot(time, position_error(1:2,:)', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Error (degree)');
title('RIEKF: Horizontal Error');
grid on;
figure;
plot(time, position_error(3,:)', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Altitude Error (m)');
title('RIEKF: Altitude Error');
grid on;
%%
function Df = f_jacobian(x,u)
% Jacobian of dynamics of remaining states [p;v]
phi = x(1);
h = x(3);
v1 = x(4);
v2 = x(5);
v3 = x(6);
ba1 = x(7);
ba2 = x(8);
ba3 = x(9);
Sa1 = x(10);
Sa2 = x(11);
Sa3 = x(12);
aB1 = u(4);
aB2 = u(5);
aB3 = u(6);


Df1 = [-(8750552118621148022459775211653*v1*pi*cos((pi*phi)/180)*sin((pi*phi)/180))/(12379400392853802748991242240*(h + 3401313089569521/(536870912*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2)))^2*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(5/2)), 0,                                                                                                                                                                                                     -v1/(h + 3401313089569521/(536870912*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2)))^2,                                                             1/(h + 3401313089569521/(536870912*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2))),                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                   0;
      (v2*pi*sin((pi*phi)/180))/(180*cos((phi*pi)/180)^2*(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))) - (16409021683231658168141*v2*pi*sin((pi*phi)/180))/(69175290276410818560*(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))^2*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2)), 0,                                                                                                                                                          -v2/(cos((pi*phi)/180)*(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))^2),                                                                                                                                                                       0,                                                                                                                                                                                                1/(cos((pi*phi)/180)*(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))),                                                                                                                                                   0;
     0, 0,                             0,                          0,            0,     -1;
                                                       v2*sin((pi*phi)/180)*((v2*pi*sin((pi*phi)/180))/(180*(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))) + (16409021683231658168141*v2*pi*cos((pi*phi)/180)^2*sin((pi*phi)/180))/(69175290276410818560*(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))^2*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2))) - (v2*pi*cos((pi*phi)/180)*((v2*cos((pi*phi)/180))/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2)) + 5380631756738151/36893488147419103232))/180 - (8750552118621148022459775211653*v1*v3*pi*cos((pi*phi)/180)*sin((pi*phi)/180))/(12379400392853802748991242240*(h + 3401313089569521/(536870912*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2)))^2*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(5/2)), 0,                                                                                                                        (v2^2*cos((pi*phi)/180)*sin((pi*phi)/180))/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))^2 - (v1*v3)/(h + 3401313089569521/(536870912*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2)))^2,                                                            v3/(h + 3401313089569521/(536870912*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2))), - sin((pi*phi)/180)*((v2*cos((pi*phi)/180))/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2)) + 5380631756738151/36893488147419103232) - (v2*cos((pi*phi)/180)*sin((pi*phi)/180))/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2)),                                        v1/(h + 3401313089569521/(536870912*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2)));
                         (v1*pi*cos((pi*phi)/180)*((v2*cos((pi*phi)/180))/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2)) + 5380631756738151/36893488147419103232))/180 - (1793543918912717*v3*pi*sin((pi*phi)/180))/2213609288845146193920 - v1*sin((pi*phi)/180)*((v2*pi*sin((pi*phi)/180))/(180*(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))) + (16409021683231658168141*v2*pi*cos((pi*phi)/180)^2*sin((pi*phi)/180))/(69175290276410818560*(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))^2*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2))) - (16409021683231658168141*v2*v3*pi*cos((pi*phi)/180)*sin((pi*phi)/180))/(69175290276410818560*(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))^2*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2)), 0,                                                                                                                                          - (v2*v3)/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))^2 - (v1*v2*cos((pi*phi)/180)*sin((pi*phi)/180))/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))^2, sin((pi*phi)/180)*((v2*cos((pi*phi)/180))/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2)) + 5380631756738151/36893488147419103232),                                                                                   v3/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2)) + (v1*cos((pi*phi)/180)*sin((pi*phi)/180))/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2)), (5380631756738151*cos((pi*phi)/180))/36893488147419103232 + v2/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2));
(1793543918912717*v2*pi*sin((pi*phi)/180))/2213609288845146193920 - (2513388269117674108869746414689*pi*cos((pi*phi)/90)*sin((pi*phi)/90))/1993841993677373809355710590420516864 + (6731709636737046437443812992459*pi*cos((pi*phi)/180)*sin((pi*phi)/180))/11682667931703362164193616740745216 + (5319273606304369*h*pi*cos((pi*phi)/180)*sin((pi*phi)/180))/108803323765316625723555840 + (16409021683231658168141*v2^2*pi*cos((pi*phi)/180)*sin((pi*phi)/180))/(69175290276410818560*(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))^2*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2)) + (8750552118621148022459775211653*v1^2*pi*cos((pi*phi)/180)*sin((pi*phi)/180))/(12379400392853802748991242240*(h + 3401313089569521/(536870912*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2)))^2*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(5/2)), 0, (713053462628379*h)/4951760157141521099596496896 + v2^2/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2))^2 + v1^2/(h + 3401313089569521/(536870912*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2)))^2 + (5319273606304369*sin((pi*phi)/180)^2)/1208925819614629174706176 - 7290625494578301/2361183241434822606848,                                                       -(2*v1)/(h + 3401313089569521/(536870912*(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(3/2))),                                                                                                                                                 - (5380631756738151*cos((pi*phi)/180))/36893488147419103232 - (2*v2)/(h + 6378137/(1 - (7718094648906879*sin((phi*pi)/180)^2)/1152921504606846976)^(1/2)),                                                                                                                                                   0];
 
Df2 = [ 0,       0,      0,                  0,                  0,                  0;
      0,       0,      0,                  0,                  0,                  0;
      0,       0,      0,                  0,                  0,                  0;
   -1/Sa1,     0,      0,   -(aB1 - ba1)/Sa1^2,           0,                  0;
       0,  -1/Sa2,    0,                  0,   -(aB2 - ba2)/Sa2^2,           0;
       0,       0, -1/Sa3,                 0,                  0, -(aB3 - ba3)/Sa3^2];
Df = [Df1 Df2; zeros(6,12)];
end

function Dh = h_jacobian(rs)
% function of remaining state
Dh = [eye(3) zeros(3,3) zeros(3, 6)];
end

function xdot = qrv_dynamics(t,x,uB,IMU)
q = x(1:4);
bg = IMU(1:3);
ba = IMU(4:6);
Sa = IMU(7:9);
omegaBhat = uB(1:3);
omegaB = omegaBhat - bg;
aBhat = uB(4:6);
aB = (1./Sa).*(aBhat - ba);
qdot = 1/2 * qmult(q,[0;omegaB]);
pdot = compute_pdot(x);
vdot = compute_vdot(x,aB);
xdot = [qdot;pdot;vdot];
end

function [Qq,Qr] = compute_Q(x,dt)
n = length(x);
nr = n-4;
sig_v_gyro = deg2rad(12e-3)*sqrt(1/dt); % 12 mdps/sqrt(Hz)
sig_u_gyro = deg2rad(3/3600); 
Qq = [sig_v_gyro*eye(3) zeros(3,3);
      zeros(3,3) sig_u_gyro*eye(3)]; % min error state dim
blkdiag(0.1*eye(4),0.1*eye(6),1e-6);
sig_v = dt * 0.03*sqrt(1/3600);
sig_p = sig_v*sqrt(dt^3/3);
sig_u_acc = 40e-6 *9.8;
nla = 5e-3; % 0.5% of +-16g  
Qr = blkdiag(sig_p*eye(3), sig_v*eye(3), sig_u_acc*eye(3), nla^2*eye(3) );
end

function [Rq,Rr] = get_R(x)
sig_mag = 1e-8;
Rq = sig_mag * eye(3); % min error state dim
Rr = diag([(1.35e-5/2)^2;(1.63e-5/2)^2;(2/2)^2]);
end