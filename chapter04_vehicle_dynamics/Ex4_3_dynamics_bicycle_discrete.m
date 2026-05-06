%Example 4.3 comparison between continuous and discrete time perturbation
%models

clc; clear; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

%% Vehicle Parameters
m = 1500;       % Vehicle mass (kg)
Lf = 1.2;       % Distance from CG to front axle (m)
Lr = 1.6;       % Distance from CG to rear axle (m)
Iz = m * ((0.5 * (Lf + Lr))^2);  % Yaw moment of inertia (kg*m^2)
Cd = 0.32;      % Drag coefficient
A = 2.2;        % Frontal area (m^2)
rho = 1.225;    % Air density (kg/m^3)
CxA = 0.5 * Cd * A * rho;  % Aerodynamic resistance coefficient
Cxf = 60000;  % Lower front longitudinal stiffness
Cxr = 120000; % Higher rear longitudinal stiffness for RWD
Cyf = 30000;    % Lateral stiffness front (N/rad)
Cyr = 30000;    % Lateral stiffness rear (N/rad)

%% Simulation Time
Ts = 0.05;  % Sampling time
T_final = 20;
N = floor(T_final / Ts);  % Number of steps
Tspan = (0:N-1) * Ts;  % Time steps

%% Control Inputs
acceleration = 2; % m/s^2
delta = 5 * pi / 180; % Steering angle (rad)
Sf = 0;
Sr = m * acceleration / (2 * Cxr);
U = [delta; Sf; Sr];
delta_nom = delta;
Sf_nom = Sf;
Sr_nom = Sr;
%% Initial Conditions
X0_nominal = [0; 0; 0; 15; 0; 0];
X0_linear = @(t) 1e-1 * [sin(0.1*t); cos(0.1*t); 0.1*sin(0.05*t); 1; 1; 1];
X0_perturbed = X0_nominal + X0_linear(0);

%% Continuous-Time Simulation
[t_nl, X_nl] = ode45(@(t, X) dynamic_bicycle_model(t, X, U, m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr), Tspan, X0_nominal);
[t_lin_c, X_lin_c] = ode45(@(t, X) linearized_bicycle_dynamics(t, X, U-U, m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr, ...
    interp1(t_nl, X_nl(:,3), t), ... % psi_nom
    interp1(t_nl, X_nl(:,4), t), ... % vx_nom
    interp1(t_nl, X_nl(:,5), t), ... % vy_nom
    interp1(t_nl, X_nl(:,6), t), ... % omega_nom
    delta_nom, Sf_nom), Tspan, X0_linear(0));

%% Discrete-Time Storage
X_nl_d = zeros(6, N);
X_lin_d = zeros(6, N);
X_nl_d(:,1) = X0_nominal;
X_lin_d(:,1) = X0_linear(0);

%% Discrete-Time Simulation Loop
for k = 1:N-1

    % Update Nominal Parameters
    psi_nom = X_nl_d(3,k);
    vx_nom = X_nl_d(4,k);
    vy_nom = X_nl_d(5,k);
    omega_nom = X_nl_d(6,k);

    
    % Recompute Transformation Coefficients based on updated nominal parameters
    alpha13 = -(vx_nom * sin(psi_nom) + vy_nom * cos(psi_nom));
    alpha14 = cos(psi_nom);
    alpha15 = -sin(psi_nom);
    alpha23 = vx_nom * cos(psi_nom) - vy_nom * sin(psi_nom);
    alpha24 = sin(psi_nom);
    alpha25 = cos(psi_nom);
    alpha36 = 1;

    alpha44 = (1/m) * (-2 * Cyf * sin(delta_nom) * (vy_nom + Lf * omega_nom) / vx_nom^2 - 2 * CxA * vx_nom);
    alpha45 = (1/m) * (m * omega_nom + 2 * Cyf * sin(delta_nom) / vx_nom);
    alpha46 = (1/m) * (m * vy_nom - 2 * Cyf * sin(delta_nom) * Lf / vx_nom);
    beta41 = (1/m) * (-2 * Cxf * Sf_nom * sin(delta_nom) - 2 * Cyf * (sin(delta_nom) + delta_nom * cos(delta_nom)));
    beta42 = (1/m) * (2 * Cxf * cos(delta_nom));
    beta43 = (1/m) * (2 * Cxr);

    alpha54 = (1/m) * (-m * omega_nom + 2 * Cyf * cos(delta_nom) * (vy_nom + Lf * omega_nom) / vx_nom^2 + 2 * Cyr * (vy_nom - Lr * omega_nom) / vx_nom^2);
    alpha55 = (1/m) * (-2 * Cyf * cos(delta_nom) / vx_nom - 2 * Cyr / vx_nom);
    alpha56 = (1/m) * (-m * vx_nom - 2 * Cyf * cos(delta_nom) * Lf / vx_nom + 2 * Cyr * Lr / vx_nom);
    beta51 = (1/m) * (2 * Cyf * (cos(delta_nom) - delta_nom * sin(delta_nom)) + 2 * Cxf * Sf_nom * cos(delta_nom));
    beta52 = (1/m) * (2 * Cxf * sin(delta_nom));

    alpha64 = (1/Iz) * (2 * Cyf * Lf * cos(delta_nom) * (vy_nom + Lf * omega_nom) / vx_nom^2 - 2 * Cyr * Lr * (vy_nom - Lr * omega_nom) / vx_nom^2);
    alpha65 = (1/Iz) * (-2 * Cyf * Lf * cos(delta_nom) / vx_nom + 2 * Cyr * Lr / vx_nom);
    alpha66 = (1/Iz) * (-2 * Cyf * Lf^2 * cos(delta_nom) / vx_nom - 2 * Cyr * Lr^2 / vx_nom);
    beta61 = (1/Iz) * (2 * Cyf * Lf * (cos(delta_nom) - delta_nom * sin(delta_nom)) + 2 * Cxf * Sf_nom * Lf * cos(delta_nom));
    beta62 = (1/Iz) * (2 * Cxf * Lf * sin(delta_nom));

    % Recompute Ad and Bd
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

    Ad = eye(size(A)) + Ts * A;
    Bd = Ts * B;

    X_nl_d(:,k+1) = X_nl_d(:,k) + Ts * dynamic_bicycle_model(0, X_nl_d(:,k), U, m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr);
    X_lin_d(:,k+1) = Ad * X_lin_d(:,k) + Bd * (U - [delta_nom; Sf_nom; Sr_nom]);
end



%% Plot Perturbation Model Comparison
figure;
subplot(3,1,1); hold on;
plot(t_lin_c, X_lin_c(:,1), 'r-', 'LineWidth', 2);
plot(Tspan, X_lin_d(1,:), 'b--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\tilde{x}$ (m)', 'Interpreter', 'latex');
legend('Continuous-Time Linearized', 'Discrete-Time Linearized');
grid on;box on;

subplot(3,1,2); hold on;
plot(t_lin_c, X_lin_c(:,2), 'r-', 'LineWidth', 2);
plot(Tspan, X_lin_d(2,:), 'b--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\tilde{y}$ (m)', 'Interpreter', 'latex');
legend('Continuous-Time Linearized', 'Discrete-Time Linearized');
grid on;box on;

subplot(3,1,3); hold on;
plot(t_lin_c, X_lin_c(:,3) * (180/pi), 'r-', 'LineWidth', 2);
plot(Tspan, X_lin_d(3,:) * (180/pi), 'b--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\tilde{\psi}$ (deg)', 'Interpreter', 'latex');
legend('Continuous-Time Linearized', 'Discrete-Time Linearized');
grid on;box on;
