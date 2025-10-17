function [x_plus,P_plus,Pq_plus] = update_mag(x_minus,P_minus,Pq_minus,Hq,Rq,R,magI,mag_meas)
    % STEP 2: KALMAN GAIN
    Kq = Pq_minus * Hq' / (Hq * Pq_minus * Hq' + Rq);
    H = [skew(magI) zeros(3,18)];
    K = P_minus*H'/(H*P_minus*H' + R);

    % STEP 3: MEAS. UPDATE
    q_minus = x_minus(1:4);
    dv = magI - q2v(sandwich(q_minus,[0;mag_meas]));   
    err = Kq*dv;
    cq = err(1:3);
    cbg = err(4:6);
    dq = quaternion_exp(-cq/2);
    % q_plus = delta_q ⊗ q_minus
    q_plus = qmult(dq,q_minus);
    bg_plus = x_minus(11:13) - q2v(sandwich(qconj(q_plus),[0;cbg]));
    
    x_plus = x_minus;
    x_plus(1:4) = q_plus;
    x_plus(11:13) = bg_plus;
    Pq_plus = (eye(6) - Kq * Hq) * Pq_minus * (eye(6) - Kq * Hq)' + Kq * Rq * Kq'; 
    P_plus = (eye(21,21) - K*H)*P_minus* (eye(21,21)-K * H)' + K*R*K';
end