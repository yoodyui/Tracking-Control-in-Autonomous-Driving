function dXdt = kinematic_bicycle_model(~, X, U, L, vx)
    % Extract States
    x = X(1);
    y = X(2);
    psi = X(3);

    % Extract Control Input
    delta = U(1); % Steering angle

    % Kinematic Bicycle Model Equations
    dxdt = vx * cos(psi);
    dydt = vx * sin(psi);
    dpsidt = (vx / L) * tan(delta);
    
    % Return State Derivative
    dXdt = [dxdt; dydt; dpsidt];
end
