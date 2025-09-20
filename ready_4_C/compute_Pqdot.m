function Pqdot = compute_Pqdot(x,Pq,Qq,w_meas)
% RIEKF error space covariance derivative
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