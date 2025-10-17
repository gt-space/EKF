
% DEFINE CONSTANTS %
we = 7.29211*10^(-5); % Earth's sidereal rotation rate [rad/s]
% DEFINE PARAMETERS
IMU_freq = 1000; % save sample rate variable
dt = 1/IMU_freq; % get time step

R = diag([1.35e-5;1.65e-5;2]); % MODIFY TO MODEL GNSS NOISE
Rq = diag([3.2e-7 4.1e-7 3.2e-7]);
Rb = 2.5e-3; % barometer noise variance


% define Q ( eqn 7.76b)
% MODIFY TO MODEL GYRO/ACCEL NOISE
nu_gv_mat = 12e-3*eye(3); % using rad/s
nu_gu_mat = 3/3600*eye(3); % using rad/s
nu_av_mat = (200e-6*9.81)*eye(3); % using m/s^2
nu_au_mat = (40e-6 *9.8)*eye(3); % using m/s^2

Q = dt*[nu_gv_mat, zeros(3,9);...
    zeros(3,3), nu_gu_mat, zeros(3,6);...
    zeros(3,6), nu_av_mat, zeros(3,3);...
    zeros(3,9), nu_au_mat];
Qq = dt*[nu_gv_mat zeros(3,3);
      zeros(3,3) nu_gu_mat];

% DEFINE INITIAL CONDITIONS %
lla0 = [35.3478813171387; -117.806831359863; 625];
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

x_prev = single(x0);
P_prev = single(P0);
Pq_prev = single(Pq0);

%% EKF
a_meas = [19.7918302847682;0.0131653560858343;-0.000292017931905962];
w_meas = [0.000328697971422360;0.000155072323595796;-0.00202279818123875];
% for i=1:50
% tic
    w_hat = compute_what(x_prev,w_meas,we); % compute gyro estimate (accounting for sf & bias)
    ahat_n = compute_ahat(x_prev,a_meas); % compute accel estimate (accounting for sf, bias & frame transformation)
    [x_minus,P_minus,Pq_minus] = propagate(x_prev,P_prev,Pq_prev,w_hat,ahat_n,w_meas,a_meas,we,Q,Qq,dt);
% toc
% end

