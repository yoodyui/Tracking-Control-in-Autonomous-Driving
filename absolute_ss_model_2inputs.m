function dXdt = absolute_ss_model_2inputs(~, X, U, L)
    % Extract States
    psi0 = X(3); % Heading angle
    v0 = X(4);   % Longitudinal velocity
    
    % Extract Control Inputs
    a0 = U(1);     % Longitudinal acceleration
    delta0 = U(2); % Steering angle

    % State-Space Matrices (Nominal Model)
    A = [0, 0, -v0 * sin(psi0), cos(psi0);
         0, 0,  v0 * cos(psi0), sin(psi0);
         0, 0,  0, (1 / L) * tan(delta0);
         0, 0,  0, 0];

    B = [0, 0;
         0, 0;
         0, (v0 / L) * sec(delta0)^2;
         1, 0];

    d = [v0 * psi0 * sin(psi0);
        -v0 * psi0 * cos(psi0);
        -(v0 / L) * delta0 * sec(delta0)^2;
         0];

    % Compute state derivative
    dXdt = A * X + B * U + d;
end
