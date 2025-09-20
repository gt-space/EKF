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