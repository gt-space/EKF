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
