function [x_plus,P_plus,Pq_plus] = update_baro(x_minus,P_minus,Pq_minus,Hb_func,Rb,press_meas,hb)
    % STEP 2: KALMAN GAIN
    Hb = Hb_func(x_minus);
    K = P_minus*Hb'/(Hb*P_minus*Hb' + Rb);

    % STEP 3: MEAS. UPDATE
    q_minus = x_minus(1:4);
    Delta_x = K*(press_meas - hb(x_minus)); % Step 3.1: compute DeltaX vector
    % q_plus = delta_q ⊗ q_minus
    q_plus = q_minus;
    p_plus = x_minus(5:7) + Delta_x(4:6);
    v_plus = x_minus(8:10) + Delta_x(7:9);
    bg_plus = x_minus(11:13);
    ba_plus = x_minus(14:16) + Delta_x(13:15);
    kg_plus = x_minus(17:19) + Delta_x(16:18);
    ka_plus = x_minus(20:22) + Delta_x(19:21);

    x_plus = [q_plus; p_plus; v_plus; bg_plus; ba_plus; kg_plus; ka_plus];
    P_plus = (eye(21,21) - K*Hb)*P_minus* (eye(21,21)-K * Hb)' + K*Rb*K';
    Pq_plus = Pq_minus;
end