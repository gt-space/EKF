function [X,IMU,Pq,Pr] = hybrid_full_RIEKF_offline(f,Df,h,Dh,x0,IMU0,Pq0,Pr0,time,U,Y,...
                                      get_Q,get_R,get_inert_mag_vec,ratio)
n = length(x0);     
nr = n - 4; % remaining state
N = length(time);
X = zeros(n,N);
X(:,1) = x0;
IMU = zeros(9,N);
IMU(:,1) = IMU0;
Pq = zeros(6,6,N);
Pq(:,:,1) = Pq0;
Pr = zeros(nr+6,nr+6,N);
Pr(:,:,1) = Pr0;

for i = 1:N-1
    [X(:,i+1),IMU(:,i+1),Pq(:,:,i+1),Pr(:,:,i+1)]... 
    = hybrid_full_RIEKF_propagate(f,Df,X(:,i),IMU(:,i),Pq(:,:,i),Pr(:,:,i),...
      [time(i) time(i+1)],U(:,i),get_Q);
    % update every interval starting from first
    if mod(i-1, ratio) == 0 && ~sum(isnan(Y(:,i)))
        [X(:,i+1),IMU(:,i+1),Pq(:,:,i+1),Pr(:,:,i+1)]... 
        = hybrid_full_RIEKF_update(h,Dh,X(:,i+1),IMU(:,i+1),Pq(:,:,i+1),Pr(:,:,i+1),...
          Y(:,i),get_R,get_inert_mag_vec);
    end
end

end