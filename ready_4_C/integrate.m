function [x_minus,P_minus,Pq_minus] = integrate(x,P,Pq,dt,qdot,pdot,vdot,Pdot,Pqdot)
% Euler integration
% assemble xdot vector
xdot = [qdot; pdot; vdot; zeros(12,1)];

x_minus = x + dt*xdot; % use Euler integration
q_next = x_minus(1:4); % get propagated quaternion

q_norm = norm(q_next);
if q_norm > 1.01 % if quaternion norm is above 1% greater than unity
%     warning('Quaternion norm of %.2f is > 1 percent of unity, consider decreasing integrator step size.',norm(q_next))
end
x_minus(1:4) = q_next/q_norm; % normalize quaternion

P_minus = P + dt*Pdot; % use Euler integration
Pq_minus = Pq + dt*Pqdot; 

% % assemble xdot vector
% xdot = [pdot; vdot; zeros(12,1)];
% q = x(1:4);
% q_next = qmult(q, quaternion_exp(qdot)); % use exponential map to get unit quaternion
% q_norm = norm(q_next);
% if q_norm > 1.01 % if quaternion norm is above 1% greater than unity
% %     warning('Quaternion norm of %.2f is > 1 percent of unity, consider decreasing integrator step size.',norm(q_next))
% end
% q_minus = q_next/q_norm; % normalize quaternion
% 
% x_minus = [q_minus; x(5:end) + dt*xdot];
% P_minus = P + dt*Pdot; % use Euler integration
% Pq_minus = Pq + dt*Pqdot; % use Euler integration

end