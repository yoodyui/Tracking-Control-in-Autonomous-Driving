function dXdt = dynamic_bicycle_model(~, X, U, m, Iz, Lf, Lr, CA, Cxf, Cxr, Cyf, Cyr)
    % Extract States
    x = X(1);
    y = X(2);
    psi = X(3);
    vx = X(4);
    vy = X(5);
    omega = X(6);

    % Extract Control Inputs
    delta = U(1); % Steering angle
    Sf = U(2);    % Front slip ratio
    Sr = U(3);    % Rear slip ratio
    
    vx_min = 0.5; % (m/s), prevent division by very small numbers

    % Compute Forces
    Fyf = Cyf * (delta - (vy + Lf * omega) / max(abs(vx), vx_min)); % Lateral force front
    Fyr = Cyr * ((-vy + Lr * omega) / max(abs(vx), vx_min));       % Lateral force rear
    Fx = 2 * Cxf * Sf * cos(delta) + 2 * Cxr * Sr - 2 * Fyf * sin(delta); % Total longitudinal force

    % Compute Nonlinear Dynamics
    dxdt = vx * cos(psi) - vy * sin(psi);
    dydt = vx * sin(psi) + vy * cos(psi);
    dpsidt = omega;
    dvxdt = (1/m) * (m * vy * omega + Fx - CA * vx^2);
    dvydt = (1/m) * (-m * vx * omega + 2 * Fyf * cos(delta) + 2 * Cxf * Sf * sin(delta) + 2 * Fyr);
    domegadt = (1/Iz) * (2 * Fyf * Lf * cos(delta) + 2 * Cxf * Sf * Lf * sin(delta) - 2 * Fyr * Lr);

    % Return State Derivative
    dXdt = [dxdt; dydt; dpsidt; dvxdt; dvydt; domegadt];
end
