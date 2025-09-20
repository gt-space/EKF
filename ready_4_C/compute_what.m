function w_hat = compute_what(x,w_meas,we)

sf_g = x(17:19); % get gyro scale factor
bias_g = x(11:13); % get gyro bias
q = x(1:4); % get quaternion

D_bn = quaternion2DCM(q)'; % get DCM I2B from quaternion, inertial to body

w_n = compute_wn(x,we);

w_hat = 1./(1 + sf_g).*(w_meas - bias_g) - D_bn*w_n; % Eqn. 7.68b

end