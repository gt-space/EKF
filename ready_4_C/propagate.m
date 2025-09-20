function [x_minus,P_minus,Pq_minus] = propagate(x_plus,P_plus,Pq_plus,w_hat,ahat_n,w_meas,a_meas,we,Q,Qq,dt)
    % step 1.5: time update: propagate to next time step
    qdot = compute_qdot(x_plus,w_hat); % compute attitude derivative
    pdot = compute_lla_dot(x_plus); % compute position derivative
    vdot = compute_vdot(x_plus,ahat_n,we); % compute velocity derivative
    Pdot = compute_Pdot(x_plus,P_plus,Q,a_meas,w_meas,we); % compute covariance derivative
    Pqdot = compute_Pqdot(x_plus,Pq_plus,Qq,w_meas);
    [x_minus,P_minus,Pq_minus] = integrate(x_plus,P_plus,Pq_plus,dt,qdot,pdot,vdot,Pdot,Pqdot);
end