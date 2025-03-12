clear
syms q [4 1] real
syms phi lambda h real
syms v [3 1] real
syms t real
syms aB [3 1] real
syms omegaB [3 1] real
syms ba [3 1] real
syms Sa [3 1] real

x = [phi;lambda;h;v];
uB = [omegaB;aB];
IMU = [ba;Sa];
Df = jacobian([qrv_dynamics(t,x,q,uB,IMU)],[x]) 
Df2 = jacobian([qrv_dynamics(t,x,q,uB,IMU)],[IMU]) 


function xdot = qrv_dynamics(t,x,q,uB,IMU)
aBhat = uB(4:6);
ba = IMU(1:3);
Sa = IMU(4:6);
aB = (1./Sa).*(aBhat - ba);
pdot = compute_pdot([q;x]);
vdot = compute_vdot([q;x],aB);
xdot = [pdot;vdot];
end
