function [xhat_plus,Phat_plus] = hybrid_EKF_update(f,Df,h,Dh,xhat,Phat,tspan,y,Qtil,Rtil)
% f: dynamics, Df: Jacobian, h: measurement, xhat: previous state estimate, 
% Phat: previous covariance estimate, Qtil: discretized process noise, 
% Rtil: discretized sensor noise,
% get_yhat: get current sensor input

n = size(xhat,1);
I = eye(n);
A = Df(xhat);
% state propagation
[~,X] = ode45(@(t,x)f(t,x),tspan,xhat);
xhat_minus = X(end,:)'; 
% covariance propagation P_{k+1}^-
phat = Phat(:);
[~,P] = ode45(@(t,P)RiccatiHybrid(t,P,A,Qtil),tspan,phat);
P_minus = reshape(P(end,:), size(Qtil));
% compute Kalman gain
H = Dh(xhat_minus);
K = P_minus*H'/(H*P_minus*H'+Rtil);
% measurement update
dy = y-h(xhat_minus);
% new state estimate
xhat_plus = xhat_minus + K*dy;
q_plus = xhat_plus(1:4);
xhat_plus(1:4) = q_plus/norm(q_plus);
% new covariance estimate
Phat_plus = (I-K*H)*P_minus*(I-K*H)'+K*Rtil*K';

function dPdt = RiccatiHybrid(t,P,A,Qtil)
    P = reshape(P,size(Qtil));
    dPdt = A*P + P*A' + Qtil;
    dPdt = dPdt(:);
end


function y = correct_IMU(S,b,yhat)
% S: IMU scaling factor, a vector representing diagonal matrix 
% b: bias
% yhat = Sy + b + v
    y = (1./S).*yhat - b;
end
end