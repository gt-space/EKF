f = @(x) 0;
Df = @(x) 0;
h = @(x) 0;
Dh = @(x) 0;
x0 = get_initial_state;
n = size(x0,1);
P0 = get_initial_covariance;
m = size(y,1);
N = 1000;
get_yhat = @()
Xhat = zeros(n,N)
while true
    [xhat,Phat] = hybrid_EKF_update(f,Df,h,Dh,xhat,Phat,Qtil,Rtil,get_yhat)
end

%%
function [xhat_plus,Phat_plus] = hybrid_EKF_update(f,Df,h,Dh, xhat,Phat,Qtil,Rtil,get_yhat)
% f: dynamics, Df: Jacobian, h: measurement, xhat: previous state estimate, 
% Phat: previous covariance estimate, Qtil: discretized process noise, 
% Rtil: discretized sensor noise,
% get_yhat: get current sensor input

n = size(xhat,1);
I = eye(n);
A = Df(xhat);
% state propagation
[~,X] = ode45(@(t,x)f(t,x),[0,dt],xhat);
xhat_minus = X(end,:)'; 
% covariance propagation P_{k+1}^-
phat = Phat(:);
[~,P] = ode45(@(t,P)RiccatiHybrid(t,P,A,Qtil),[0,dt],phat);
P_minus = reshape(P(end,:), size(Qtil));
% compute Kalman gain
H = Dh(xhat_minus);
K = P_minus*H'/(H*P_minus*H'+Rtil);
% measurement update
yhat = get_yhat;
dy = yhat-h(xhat_minus);
% new state estimate
xhat_plus = xhat_minus + K*dy;
% new covariance estimate
Phat_plus = (I-K*H)*P_minus*(I-K*H)'+K*Rtil*K';
end

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