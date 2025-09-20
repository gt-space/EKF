function ahat_n = compute_ahat(x,a_meas)

sf_a = x(20:22); % get accel scale factor
bias_a = x(14:16); % get accel bias

ahat_b = 1./(1 + sf_a).*(a_meas - bias_a); % accel. estimate in body frame (eqn. 7.68j)

D_nb = quaternion2DCM(x(1:4)); % get DCM from body frame to NED frame 

ahat_n = D_nb*ahat_b; % get accel. estimate in NED frame
end