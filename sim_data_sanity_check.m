clear;
load("raw_flight_sim_data.mat");
% apogee_index = 9447;
final_index = 104932; % 100 secs, covering apogee
% final_index = 50000; % 100 secs, covering apogee
posI = double(simData.pos(:,1:final_index)); % NED
velI = double(simData.vel(:,1:final_index)); % NED
accI = double(simData.acc(:,1:final_index)); % NED
q = double(simData.quat(:,1:final_index));
omegaB = double(simData.omega(:,1:final_index)); % body
time = double(simData.time(:,1:final_index));
dt = time(2) - time(1);
n=6;
Xhistory = zeros(n,final_index);
Xhistory(:,1) = [posI(:,1);velI(:,1)];
Yhistory(:,1) = [posI(:,1);velI(:,1)];

for i = 1:final_index-1
    acc_interp = linear_interp(accI(:,i),accI(:,i+1));
    [~,X] = ode45(@(t,x)dynamics(t,x,acc_interp),[0,dt],Xhistory(:,i));
    Xhistory(:,i+1) = X(end,:)';
    Yhistory(:,i+1) = Yhistory(:,i)+dt*dynamics(0,Yhistory(:,i),acc_interp);
end
figure
plot(time,Xhistory(1:3,:)-posI);
figure
plot(time,Yhistory(1:3,:)-posI);
figure
plot(time,Xhistory(4:6,:)-velI);
figure
plot(time,Yhistory(4:6,:)-velI);
function xdot = dynamics(t,x,acc)
pdot = x(4:6);
vdot = acc(t);
xdot = [pdot;vdot];
end
function a_interp = linear_interp(a0,a1)
a_interp = @(t) (1-t)*a0 + t*(a1);
end