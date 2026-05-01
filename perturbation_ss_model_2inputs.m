function dXdt = perturbation_ss_model_2inputs(~, X, U, L, v0, psi0, delta0)

    % State-Space Matrices for Perturbation Model (including longitudinal dynamics)
    A = [0, 0, -v0 * sin(psi0), cos(psi0);
         0, 0,  v0 * cos(psi0), sin(psi0);
         0, 0,  0, (1 / L) * tan(delta0);
         0, 0,  0, 0];

    B = [0, 0;
         0, 0;
         0, (v0 / L) * sec(delta0)^2;
         1, 0];

    % Compute state derivative (no offset term in perturbation model)
    dXdt = A * X + B * U;
end
