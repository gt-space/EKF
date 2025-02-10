function [Xhat,Phat] = hybrid_RIEKF_offline(f,Df,h,Dh,x0,P0,time,Y,get_Qtil,get_Rtil)
n = size(x0,1);
m = size(Y,1);
N = length(time);
Xhat = zeros(n,N);
Xhat(:,1) = x0;
Phat = zeros(n-1,n-1,N);
Phat(:,:,1) = P0;
%% --- Build the Transformation Matrix F_mat ---
% For a right-invariant EKF the minimal attitude error (delta_theta, 3x1) is
% mapped to an additive quaternion error by:
%    delta_q ≈ [1; 0.5*delta_theta]

T = [ zeros(1,3); 0.5 * eye(3) ];  % 4x3 mapping matrix for quaternion part
% For the remaining states, the mapping is identity.
F_mat = [ T,               zeros(4, n-4);
          zeros(n-4, 3), eye(n-4) ];  % size: n_full x (n_full-1)
%%
for i = 1:N-1
    Qtil = get_Qtil(Xhat(:,i));
    Rtil = get_Rtil(Xhat(:,i));
    [Xhat(:,i+1),Phat(:,:,i+1)] = hybrid_RIEKF_update(f,Df,h,Dh,Xhat(:,i),Phat(:,:,i),[time(i) time(i+1)],Y(:,i+1),Qtil,Rtil,F_mat);
end

%%
