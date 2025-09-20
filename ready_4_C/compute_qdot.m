function qdot = compute_qdot(x,what)
q_bn = x(1:4); % get quaternion
w_quat = [0; what]; % convert angular velocity into quaternion form
qdot = 1/2*qmult(q_bn,w_quat); % compute quaternion derivative

% % instead of computing qdot, we compute q^{-1}qdot for more accurate
% % integration using quaternion exponential map later
% qdot = 0.5 * what;
end