function dXdt = kinematic_bicycle_model_2inputs(~, X, U, L)
    % Extract States
    psi = X(3); % Heading angle
    v = X(4);   % Longitudinal velocity

    % Extract Control Inputs
    a = U(1);     % Longitudinal acceleration
    delta = U(2); % Steering angle

    % Compute State Derivatives using the Kinematic Bicycle Model
    dXdt = [v * cos(psi);         % dx/dt
            v * sin(psi);         % dy/dt
            (v / L) * tan(delta); % dψ/dt
            a];                   % dv/dt
end
