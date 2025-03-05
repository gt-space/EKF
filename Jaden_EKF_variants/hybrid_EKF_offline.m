function [Xhat,Phat] = hybrid_EKF_offline(f,Df,h,Dh,x0,P0,time,Y,get_Qtil,get_Rtil)
n = size(x0,1);
m = size(Y,1);
N = length(time);
Xhat = zeros(n,N);
Xhat(:,1) = x0;
Phat = zeros(n,n,N);
Phat(:,:,1) = P0;
for i = 1:N-1
    Qtil = get_Qtil(Xhat(:,i));
    Rtil = get_Rtil(Xhat(:,i));
    [Xhat(:,i+1),Phat(:,:,i+1)] = hybrid_EKF_update(f,Df,h,Dh,Xhat(:,i),Phat(:,:,i),[time(i) time(i+1)],Y(:,i+1),Qtil,Rtil);
end

%%
