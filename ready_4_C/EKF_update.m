function [x_plus,P_plus,Pq_plus,deploy_drogue,deploy_chute] = EKF_update(x_prev,P_prev,Pq_prev,Q,Qq,H,Hq,hb,Hb,R,Rq,Rb,a_meas,w_meas,lla_meas,mag_meas,press_meas,we,magI,dt)

    w_hat = compute_what(x_prev,w_meas,we); % compute gyro estimate (accounting for sf & bias)
    ahat_n = compute_ahat(x_prev,a_meas); % compute accel estimate (accounting for sf, bias & frame transformation)

    [x_minus,P_minus,Pq_minus] = propagate(x_prev,P_prev,Pq_prev,w_hat,ahat_n,w_meas,a_meas,we,Q,Qq,dt);

    x_plus = x_minus;
    P_plus = P_minus;
    Pq_plus = Pq_minus;
    if ~sum(isnan(lla_meas))
        [x_plus,P_plus,Pq_plus] = update_GPS(x_plus,P_plus,Pq_plus,H,R,lla_meas);
    end
    if ~sum(isnan(mag_meas))
        [x_plus,P_plus,Pq_plus] = update_mag(x_plus,P_plus,Pq_plus,Hq,Rq,R,magI,mag_meas);
    end
    if ~sum(isnan(press_meas))
        [x_plus,P_plus,Pq_plus] = update_baro(x_plus,P_plus,Pq_plus,Hb,Rb,press_meas,hb);
    end

    % velocity in NED, down is positive; position in lla, down is negative
    if x_plus(10) > 0 && x_prev(7) > x_plus(7)
        deploy_drogue = 1;
    else
        deploy_drogue = 0;
    end
    % TODO: integrate terrain map
    deploy_chute = 0;
    if any(diag(P_plus)<0)
        P_plus = nearestPSD(P_plus); %+ diag(max(1e-6 - diag(P_plus), 0)); % heuristic to correct diagonal negative covariance
    end

end