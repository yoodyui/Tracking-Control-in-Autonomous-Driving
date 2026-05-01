%Example 4.2 comparison of nonlinear dynamics and dynamic perturbation
%model, perturbing the initial states
clc; clear; close all;

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

%% Control Inputs (Constant Inputs for Demonstration)
a = 2;                % Acceleration (m/s^2)
delta = 5 * pi / 180; % Steering angle (rad)
Sf = 0;  % No front slip for RWD
Sr = m * a / (2 * Cxr); % Rear wheels provide acceleration
U = [delta; Sf; Sr];   % Control input vector

%% Apply Perturbation at Initial conditions
X0_nominal = [0; 0; 0; 15; 0; 0]; % Nominal initial trajectory
X0_linear = @(t) 1e-1 * [sin(0.1*t); cos(0.1*t); 0.1*sin(0.05*t); 1; 1; 1]; % Time-varying perturbation
X0_perturbed = X0_nominal + X0_linear(0); % Initial perturbed state

%% Solve Nonlinear Dynamic Model using ODE45
[t_nl, X_nl_nominal] = ode45(@(t, X) dynamic_bicycle_model(t, X, U, m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr), Tspan, X0_nominal);
[t_nl, X_nl_perturbed] = ode45(@(t, X) dynamic_bicycle_model(t, X, U, m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr), Tspan, X0_perturbed);
% Compute Perturbation (Difference)
X_perturbation = X_nl_perturbed - X_nl_nominal;

%% Solve Linearized Dynamic Model using ODE45
[t_lin, X_lin] = ode45(@(t, X) linearized_bicycle_dynamics(t, X, U-U, m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr, ...
    interp1(t_nl, X_nl_nominal(:,3), t), ... % psi_nom
    interp1(t_nl, X_nl_nominal(:,4), t), ... % vx_nom
    interp1(t_nl, X_nl_nominal(:,5), t), ... % vy_nom
    interp1(t_nl, X_nl_nominal(:,6), t), ... % omega_nom
    delta, Sf), Tspan, X0_linear(0)); % Initial linear perturbation
%% Extract Linearized Model Results
vx_lin = X_lin(:,4);
vy_lin = X_lin(:,5);
omega_lin = X_lin(:,6);

%%Plot
%% Extract Paths (X-Y Positions)
x_nl_nominal = X_nl_nominal(:,1); 
y_nl_nominal = X_nl_nominal(:,2);
x_nl_perturbed = X_nl_perturbed(:,1); 
y_nl_perturbed = X_nl_perturbed(:,2);

% Linearized perturbation path (perturbation is around the nominal trajectory)
x_lin = X_nl_nominal(:,1) + X_lin(:,1); % Nominal x + Perturbation x
y_lin = X_nl_nominal(:,2) + X_lin(:,2); % Nominal y + Perturbation y

%% Plot X-Y Path Comparison
figure;hold on;
plot(x_nl_perturbed, y_nl_perturbed, 'r:', 'LineWidth', 2);
plot(x_lin, y_lin, 'b--', 'LineWidth', 2);
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Perturbed Nonlinear Path', 'Linearized Model');
grid on;box on;
axis equal; % Keep aspect ratio to compare paths accurately

%% Plot Perturbations in X, Y, and Psi
figure;

subplot(3,1,1); hold on;
plot(t_nl, X_perturbation(:,1), 'r:', 'LineWidth', 2);
plot(t_lin, X_lin(:,1), 'b--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\tilde{x}$ (m)', 'Interpreter', 'latex');
legend('Nonlinear Perturbation', 'Linearized Model');
grid on;box on;

subplot(3,1,2); hold on;
plot(t_nl, X_perturbation(:,2), 'r:', 'LineWidth', 2);
plot(t_lin, X_lin(:,2), 'b--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\tilde{y}$ (m)', 'Interpreter', 'latex');
legend('Nonlinear Perturbation', 'Linearized Model');
grid on;box on;

subplot(3,1,3); hold on;
plot(t_nl, X_perturbation(:,3) * (180/pi), 'r:', 'LineWidth', 2);
plot(t_lin, X_lin(:,3) * (180/pi), 'b--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\tilde{\psi}$ (deg)', 'Interpreter', 'latex');
legend('Nonlinear Perturbation', 'Linearized Model');
grid on;box on;

figure;
subplot(3,1,1); hold on;
plot(t_nl, X_perturbation(:,4), 'r:', 'LineWidth', 2);
plot(t_lin, X_lin(:,4), 'b--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\tilde{v_x}$ (m/s)', 'Interpreter', 'latex');
legend('Nonlinear Perturbation', 'Linearized Model');
grid on;box on;

subplot(3,1,2); hold on;
plot(t_nl, X_perturbation(:,5), 'r:', 'LineWidth', 2);
plot(t_lin, X_lin(:,5), 'b--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\tilde{v_y}$ (m/s)', 'Interpreter', 'latex');
legend('Nonlinear Perturbation', 'Linearized Model');
grid on;box on;

subplot(3,1,3); hold on;
plot(t_nl, X_perturbation(:,6), 'r:', 'LineWidth', 2);
plot(t_lin, X_lin(:,6), 'b--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\tilde{\omega}$ (rad/s)', 'Interpreter', 'latex');
legend('Nonlinear Perturbation', 'Linearized Model');
grid on;box on;