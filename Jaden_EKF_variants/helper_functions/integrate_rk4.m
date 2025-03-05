function [t_array, x_array] = integrate_rk4(f, t0, tf, x0, N)
% integrate_rk4 integrates the ODE dx/dt = f(t,x) from t0 to tf using RK4.
%
% Syntax:
%   [t_array, x_array] = integrate_rk4(f, t0, tf, x0, dt)
%
% Inputs:
%   f   - function handle for the ODE, f(t,x)
%   t0  - initial time
%   tf  - final time
%   x0  - initial state (column vector)
%   dt  - fixed time step
%
% Outputs:
%   t_array - vector of time points from t0 to tf
%   x_array - matrix of state values, each column corresponds to the state at a time in t_array
%
% Example:
%   % Integrate the simple ODE: dx/dt = -x, x(0)=1, dt=0.1, from t=0 to t=5
%   f = @(t,x) -x;
%   [t, x] = integrate_rk4(f, 0, 5, 1, 0.1);
%   plot(t, x);

    % Create time vector
    dt = (tf-t0)/(N-1);
    t_array = t0:dt:tf;
    
    % Preallocate state array
    % Assume state is a column vector; x_array will have as many columns as time steps.
    x_array = zeros(length(x0), N);
    x_array(:,1) = x0;
    
    % Loop over time steps
    for k = 1:N-1
        t_current = t_array(k);
        x_current = x_array(:,k);
        % RK4 step:
        x_next = rk4(f, t_current, x_current, dt);
        x_array(:, k+1) = x_next;
    end
end

