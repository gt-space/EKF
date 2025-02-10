uB = [0;0;1];
magI = [20;5;35];
J = ones(3,1);
dt = 0.01;
t0 = 0;
tf = 1;
tspan = t0:dt:tf;

h = @(x) quaternion2DCM(x(1:4))'*magI; % magnetometer measures magB
f = @(t,x) attitude_dynamics(t,x,uB,J);
Dh = @(x) [h_jacobian(x(1:4),magI), zeros(3,3)];
Df = @(x) f_jacobian(x,J);
q0 = [0;0;0;1];
omegaB0 = [0;0;0];
x0 = [q0;omegaB0];
n = size(x0,1);
P0 = 0.01*eye(n-1); % min error state dim
Q = 0.01*eye(n-1); % min error state dim
L = eye(n-1);
Qtil = L*Q*L';
R = 1*eye(3);
M = eye(3);
Rtil = M*R*M';
[time,true_state,measurements] = simulate(f,x0,tspan,h,R);
get_Qtil = @(x) Qtil;
get_Rtil = @(x) Rtil;

[Xhat,Phat] = hybrid_RIEKF_offline(f,Df,h,Dh,x0,P0,time,measurements,get_Qtil,get_Rtil);

% t= 0;
% syms q [4 1] real
% syms omegaB [3 1] real
% syms uB [3 1] real
% syms J [3 1] real
% x = [q;omegaB];
% Df = jacobian([attitude_dynamics(t,x,uB,J)],[x]);

%%
function Df = f_jacobian(x,J)
tmp = num2cell(x);
[q1, q2, q3, q4, omegaB1, omegaB2, omegaB3] = tmp{:};
tmp = num2cell(J);
[J1, J2, J3] = tmp{:};
Df = [
        0, -omegaB1/2, -omegaB2/2, -omegaB3/2,       -q2/2,       -q3/2,       -q4/2;
 omegaB1/2,         0,  omegaB3/2, -omegaB2/2,        q1/2,       -q4/2,        q3/2;
 omegaB2/2, -omegaB3/2,         0,  omegaB1/2,        q4/2,        q1/2,       -q2/2;
 omegaB3/2,  omegaB2/2, -omegaB1/2,         0,       -q3/2,        q2/2,        q1/2;
        0,         0,         0,         0,          0, (J2*omegaB3 - J3*omegaB3)/J1,  (J2*omegaB2 - J3*omegaB2)/J1;
        0,         0,         0,         0, -(J1*omegaB3 - J3*omegaB3)/J2,        0, -(J1*omegaB1 - J3*omegaB1)/J2;
        0,         0,         0,         0,  (J1*omegaB2 - J2*omegaB2)/J3,  (J1*omegaB1 - J2*omegaB1)/J3,         0
];
end
function Dh = h_jacobian(q,magI)
tmp = num2cell(q);
[q1, q2, q3, q4] = tmp{:};
tmp = num2cell(magI);
[magI1, magI2, magI3] = tmp{:};
Dh = [
    2*magI2*q4 - 2*magI3*q3,              2*magI2*q3 + 2*magI3*q4,  2*magI2*q2 - 4*magI1*q3 - 2*magI3*q1,  2*magI2*q1 - 4*magI1*q4 + 2*magI3*q2;
    2*magI3*q2 - 2*magI1*q4,  2*magI1*q3 - 4*magI2*q2 + 2*magI3*q1,              2*magI1*q2 + 2*magI3*q4,  2*magI3*q3 - 4*magI2*q4 - 2*magI1*q1;
    2*magI1*q3 - 2*magI2*q2,  2*magI1*q4 - 2*magI2*q1 - 4*magI3*q2,  2*magI1*q1 + 2*magI2*q4 - 4*magI3*q3,              2*magI1*q2 + 2*magI2*q3
];
end
function xdot = attitude_dynamics(t,x,uB,J)
% J moment of inertia as a vector
q = x(1:4);
omegaB = x(5:7);
qdot = 1/2 * qmult(q,[0;omegaB]);
tauB = uB;
% uB = sandwich(qconj(q),u);
omegaBdot = (1./J).*(tauB - cross(omegaB,J.*omegaB));
xdot = [qdot;omegaBdot];
end

function Qtil = compute_Qtil(x,Q)
% Extract quaternion components
q = x(1:4);
s  = q(1);  % Scalar part
v1 = q(2);  % First vector component
v2 = q(3);  % Second vector component
v3 = q(4);  % Third vector component

% Compute the Jacobian matrix L = partial f/partial w
L = 0.5 * [ s, -v1, -v2, -v3;
            v1,  s,  v3, -v2;
            v2, -v3,  s,  v1;
            v3,  v2, -v1,  s ];
Qtil = Q;
Qtil(1:4,1:4) = L*Q(1:4,1:4)*L';
end

function S = skew(v)
v1 = v(1);
v2 = v(2);
v3 = v(3);
S = [0 -v3 v2;
     v3 0 -v1;
     -v2 v1 0];
end

function CB2I = quaternion2DCM(q)
% quaternion describes rotation from inertial frame to body frame
% so DCM should convert body coordinates to inertial coordinates
% since frame and coordinates vary contravariantly
s = q(1);
v = q(2:4);
CB2I = (1-2*dot(v,v))*eye(3) + 2*(v*v') + 2*s*skew(v);
end

function qbar = qconj(q)
qbar = [-q(1); q(2:4)];
end

function y = sandwich(q,x)
% only for unit q
y = qmult(q,qmul(x,qconj(q)));
end

function q_out = qmult(q1, q2)
    % Quaternion multiplication: q_out = q1 * q2
    % Assumes quaternions are in [scalar, vector] form: [q0, q1, q2, q3]
    
    % Extract components
    s1 = q1(1); v1 = q1(2:4);  % Scalar and vector part of q1
    s2 = q2(1); v2 = q2(2:4);  % Scalar and vector part of q2
    
    % Compute quaternion product
    s_out = s1 * s2 - dot(v1, v2);  % Scalar part
    v_out = s1 * v2 + s2 * v1 + cross(v1, v2);  % Vector part
    
    % Combine into output quaternion
    q_out = [s_out; v_out];
end