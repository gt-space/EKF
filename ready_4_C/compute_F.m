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