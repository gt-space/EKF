function [x_plus,P_plus,Pq_plus] = update_GPS(x_minus,P_minus,Pq_minus,H,R,lla_meas)
    % STEP 2: KALMAN GAIN
    % underweighting
    if trace(P_minus(4:5,4:5)) + P_minus(6,6) > 1000
        beta = 0.25;
    else
        beta = 0;
    end
    W = (1+beta)*H*P_minus*H' + R;
    K = P_minus*H'/ W;
    % K = P_minus*H'/(H*P_minus*H' + R);

    % STEP 3: MEAS. UPDATE
    Delta_x = K*(lla_meas - x_minus(5:7)); % Step 3.1: compute DeltaX vector

    q_plus = x_minus(1:4);
    p_plus = x_minus(5:7) + Delta_x(4:6);
    v_plus = x_minus(8:10) + Delta_x(7:9);
    bg_plus = x_minus(11:13);
    ba_plus = x_minus(14:16) + Delta_x(13:15);
    kg_plus = x_minus(17:19) + Delta_x(16:18);
    ka_plus = x_minus(20:22) + Delta_x(19:21);
    x_plus = [q_plus; p_plus; v_plus; bg_plus; ba_plus; kg_plus; ka_plus];
    P_plus = (eye(21,21) - K*H)*P_minus* (eye(21,21)-K * H)' + K*R*K';
    Pq_plus = Pq_minus;   
end