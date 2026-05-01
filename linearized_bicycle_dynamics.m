function dX = linearized_bicycle_dynamics(~, X, U, m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr, psi0, vx0, vy0, omega0, delta0, Sf0)

    % Transformation Coefficients
    alpha13 = -(vx0 * sin(psi0) + vy0 * cos(psi0));
    alpha14 = cos(psi0);
    alpha15 = -sin(psi0);
    alpha23 = vx0 * cos(psi0) - vy0 * sin(psi0);
    alpha24 = sin(psi0);
    alpha25 = cos(psi0);
    alpha36 = 1;

   % Longitudinal Velocity Equations
    alpha44 = (1/m) * (-2 * Cyf * sin(delta0) * (vy0 + Lf * omega0) / vx0^2 - 2 * CxA * vx0);
    alpha45 = (1/m) * (m * omega0 + 2 * Cyf * sin(delta0) / vx0);
    alpha46 = (1/m) * (m * vy0 - 2 * Cyf * sin(delta0) * Lf / vx0);
    beta41 = (1/m) * (-2 * Cxf * Sf0 * sin(delta0) - 2 * Cyf * (sin(delta0) + delta0 * cos(delta0)));
    beta42 = (1/m) * (2 * Cxf * cos(delta0));
    beta43 = (1/m) * (2 * Cxr);

    % Lateral Velocity Equations
    alpha54 = (1/m) * (-m * omega0 + 2 * Cyf * cos(delta0) * (vy0 + Lf * omega0) / vx0^2 + 2 * Cyr * (vy0 - Lr * omega0) / vx0^2);
    alpha55 = (1/m) * (-2 * Cyf * cos(delta0) / vx0 - 2 * Cyr / vx0);
    alpha56 = (1/m) * (-m * vx0 - 2 * Cyf * cos(delta0) * Lf / vx0 + 2 * Cyr * Lr / vx0);
    beta51 = (1/m) * (2 * Cyf * (cos(delta0) - delta0 * sin(delta0)) + 2 * Cxf * Sf0 * cos(delta0));
    beta52 = (1/m) * (2 * Cxf * sin(delta0));

    % Yaw Rate Equations
    alpha64 = (1/Iz) * (2 * Cyf * Lf * cos(delta0) * (vy0 + Lf * omega0) / vx0^2 - 2 * Cyr * Lr * (vy0 - Lr * omega0) / vx0^2);
    alpha65 = (1/Iz) * (-2 * Cyf * Lf * cos(delta0) / vx0 + 2 * Cyr * Lr / vx0);
    alpha66 = (1/Iz) * (-2 * Cyf * Lf^2 * cos(delta0) / vx0 - 2 * Cyr * Lr^2 / vx0);
    beta61 = (1/Iz) * (2 * Cyf * Lf * (cos(delta0) - delta0 * sin(delta0)) + 2 * Cxf * Sf0 * Lf * cos(delta0));
    beta62 = (1/Iz) * (2 * Cxf * Lf * sin(delta0));

    
    % State-Space Representation
    A = [0 0 alpha13 alpha14 alpha15 0;
         0 0 alpha23 alpha24 alpha25 0;
         0 0 0 0 0 alpha36;
         0 0 0 alpha44 alpha45 alpha46;
         0 0 0 alpha54 alpha55 alpha56;
         0 0 0 alpha64 alpha65 alpha66];

    B = [0 0 0;
         0 0 0;
         0 0 0;
         beta41 beta42 beta43;
         beta51 beta52 0;
         beta61 beta62 0];

    % Compute Linearized State Derivative
    dX = A * X + B * U;
end
