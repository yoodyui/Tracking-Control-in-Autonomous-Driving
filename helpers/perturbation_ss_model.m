function dXdt = perturbation_ss_model(~, X, U, L, v, psi_0, delta_0)
% Continuous-time perturbation state-space model
% Inputs:
%   X - State vector [tilde_x; tilde_y; tilde_psi]
%   U - Control input [(delta - delta_0)]
%   L - Wheelbase (m)
%   v_0 - Nominal velocity (m/s)
%   psi_0 - Nominal yaw angle (rad)
%   delta_0 - Nominal steering angle (rad)
% Output:
%   dXdt - Time derivative of the perturbation state vector

%% State-Space Matrices for Perturbation Model
A = [0, 0, -v * sin(psi_0);
     0, 0,  v * cos(psi_0);
     0, 0,  0];

B = [0;
     0;
     (v / L) * sec(delta_0)^2];

% Compute state derivative (no offset term in perturbation model)
dXdt = A * X + B * U;
end
