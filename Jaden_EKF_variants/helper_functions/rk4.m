function x_next = rk4(f, t, x, dt)
% rk4 - Fixed time interval 4th-order Runge-Kutta integrator.
%
% Syntax:  x_next = rk4(f, t, x, dt)
%
% Inputs:
%    f  - Function handle representing the ODE: f(t,x)
%    t  - Current time
%    x  - Current state (column vector)
%    dt - Time step for integration
%
% Outputs:
%    x_next - State at time t+dt
%
% Example:
%    % Define the ODE: dx/dt = -x
%    f = @(t,x) -x;
%    x0 = 1;
%    dt = 0.1;
%    t0 = 0;
%    x1 = rk4(f, t0, x0, dt);

    k1 = f(t, x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt, x + dt * k3);
    
    x_next = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end
