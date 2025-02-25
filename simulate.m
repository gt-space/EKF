function [time,true_state,measurements] = simulate(f,x0,tspan,h,R)
% create simulated true trajectory and noisy measurements
n = size(x0,1);
m = size(R,1);
[time,true_state] = ode45(@(t,x)f(t,x),tspan,x0);
true_state = true_state';
N = size(true_state,2);
measurements = zeros(m,N);
for i = 1:N
    measurements(:,i) = h(true_state(1:4,i)) + R*randn(m,1);
end

end
