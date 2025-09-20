function [x_plus,P_plus,Pq_plus] = update_mag(x_minus,P_minus,Pq_minus,Hq,Rq,magI,mag_meas)
    % STEP 2: KALMAN GAIN
    Kq = Pq_minus * Hq' / (Hq * Pq_minus * Hq' + Rq);

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
    % if norm(bg_plus) > 1e-3
    %    bg_plus = zeros(3,1);
    % end
    x_plus = x_minus;
    x_plus(1:4) = q_plus;
    x_plus(11:13) = bg_plus;
    P_plus = P_minus;
    Pq_plus = (eye(6) - Kq * Hq) * Pq_minus * (eye(6) - Kq * Hq)' + Kq * Rq * Kq';   
end