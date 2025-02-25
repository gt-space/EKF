function [x_plus, IMU_plus, Pq_plus, Pr_plus] = hybrid_RIEKF_update(f, Df,...
    h, Dh, x, IMU, Pq, Pr, tspan, u, y, get_Q, get_R,get_inert_mag_vec)
    q2v = @(q) q(2:4);
    n = length(x);     
    Iq = eye(6); % error state dimension: 3 for attitude + 3 for gyro bias
    nr = n - 4;
    Ir = eye(nr+6); % including acceleration bias and scaling factor

    %% --- Nominal State Propagation ---
    % Propagate the full state using the dynamics f (e.g., with ode45)
    [~, X] = ode45(@(t,x) f(t,x), tspan, x);
    x_minus = X(end,:)';  % predicted (nominal) state
    q_minus = x_minus(1:4);
    q_minus = q_minus / norm(q_minus);  % ensure quaternion is normalized
    r_minus = x_minus(5:end);
    bg_minus = IMU(1:3); % assume zero derivative
    ba_minus = IMU(4:6); % assume zero derivative
    Sa_minus = IMU(7:9); % assume zero derivative

    %% --- Covariance Propagation (Minimal State) ---
    % Get the full state Jacobian at the previous estimate (or you could use xhat_minus)
    q = x(1:4);
    omegaB = u(1:3);
    CB2I = quaternion2DCM(q);
    I_omega = CB2I*(omegaB - bg_minus);
    F = [zeros(3,3) -eye(3);
         zeros(3,3) skew(I_omega)];
    G = [CB2I zeros(3,3);
         zeros(3,3) -CB2I];
    A = Df(r_minus,u(4:end));  % A is Jacobian of R^n part of dynamics
    pq = Pq(:);
    pr = Pr(:);
    % Integrate the Riccati equation using the minimal dynamics.
    [Qq, Qr] = get_Q(x);
    Qqtil = G*Qq*G';
    Qrtil = Qr;
    [~, Pq_sol] = ode45(@(t,P) RiccatiHybrid(t, Pq, F, Qqtil), tspan, pq);
    [~, Pr_sol] = ode45(@(t,P) RiccatiHybrid(t, Pr, A, Qrtil), tspan, pr);
    Pq_minus = reshape(Pq_sol(end,:), size(Qq));
    Pr_minus = reshape(Pr_sol(end,:), size(Qr));

    %% --- Measurement Update ---
    vI = get_inert_mag_vec(x_minus);
    Hq = [skew(vI) zeros(3,3)];
    Hr = Dh(r_minus);
    [Rqtil, Rrtil] = get_R(x_minus);
    % Compute Kalman gain in the minimal state space.
    Kq = Pq_minus * Hq' / (Hq * Pq_minus * Hq' + Rqtil);
    Kr = Pr_minus * Hr' / (Hr * Pr_minus * Hr' + Rrtil);
    
    % Compute the innovation (measurement residual)
    vB = y(1:3);
    yr = y(4:end);
    dv = vI - q2v(sandwich(q_minus,[0;vB]));
    dyr = yr - h(r_minus);
    
    err = Kq*dv;
    cq = err(1:3);
    cbg = err(4:6);
    dr = Kr * dyr;
    % Covariance update for the minimal error state.
    Pq_plus = (Iq - Kq * Hq) * Pq_minus * (Iq - Kq * Hq)' + Kq * Rqtil * Kq';
    Pr_plus = (Ir - Kr * Hr) * Pr_minus * (Ir - Kr * Hr)' + Kr * Rrtil * Kr';

    %% --- Correction Step (Quaternion Update using the Exponential Map) ---
    % Compute the correction quaternion from the minimal attitude error using the exponential map.
    dq = quaternion_exp(-cq/2);
    % q_plus = delta_q ⊗ q_minus
    q_plus = qmult(dq,q_minus);
    % q_plus = q_plus / norm(q_plus);  % re-normalize
    bg_plus = bg_minus - q2v(sandwich(qconj(q_plus),[0;cbg]));
    % Update the remaining states additively.
    r_plus = r_minus + dr(1:end-6);
    ba_plus = ba_minus + dr(end-5:end-3);
    Sa_plus = Sa_minus + dr(end-2:end);
    % Compose the full updated state.
    x_plus = [q_plus; r_plus];
    IMU_plus = [bg_plus;ba_plus;Sa_plus];
end

%% --- Helper Functions ---

function dPdt = RiccatiHybrid(t, P, A, Qtil)
    % Integrates the differential Riccati equation for the minimal error covariance.
    P = reshape(P,size(Qtil));
    dPdt = A * P + P * A' + Qtil;
    dPdt = dPdt(:);
end