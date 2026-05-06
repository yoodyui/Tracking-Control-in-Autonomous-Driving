function dXdt = absolute_ss_model(~, X, U, L, v)
% Continuous-time absolute state-space model
% Inputs:
%   X - State vector [x; y; psi]
%   A - System matrix
%   B - Input matrix
%   d - Offset term
%   delta - Steering input (constant or time-varying)
% Output:
%   dXdt - Time derivative of the state vector
psi_0 = X(3);
delta_0 = U(1); % Steering angle

%% State-Space Matrices
A = [0, 0, -v * sin(psi_0);
    0, 0,  v * cos(psi_0);
    0, 0,  0];

B = [0;
    0;
    (v / L) * sec(delta_0)^2];

d = [v * cos(psi_0) + v * psi_0 * sin(psi_0);
    v * sin(psi_0) - v * psi_0 * cos(psi_0);
    (v / L) * tan(delta_0) - (v / L) * delta_0 * sec(delta_0)^2];

% Compute state derivative
dXdt = A * X + B * U + d;
end
