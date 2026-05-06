%For Figures in Example 2, chapter 3
%The nonlinear kinematic model with 4 states and 2 inputs can be splitted into two models: absolute
%and perturbation where the input to absolute is delta_0 and that for
%perturbation is u_k-u_abs (kinematic = absolute + perturbation)
%This code needs the following function files: 
% (1) kinematic_bicycle_model_2inputs.m
% (2) absolute_ss_model_2inputs.m
% (3) perturbation_ss_model_2inputs.m
clc; clear; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

%% Vehicle Parameters
L = 3.0;         % Wheelbase (m)
a0 = 0.1;        % Nominal acceleration
psi0 = 0.1;      % Nominal yaw angle (rad)
delta0 = 5 * pi / 180;  % Nominal steering angle (rad)

%% Initial Conditions
X0_k = [0; 0; 0; 0];     % Initial state for kinematic model: [x; y; psi; v]
X0_abs = [0; 0; 0; 0];   % Initial state for absolute model: [x; y; psi; v]
X0_p = [0; 0; 0; 0];     % Initial state for perturbation model: [x; y; psi; v]

%% Simulation Time
Tspan = 0:0.1:45;  % Time vector

%% Define Input Signals
U_k = @(t) [a0 + 0.01 * sin(0.3*t); delta0 + 1*pi/180*sin(0.5*t)]; % Varying acceleration and steering
U_abs = @(t) [a0; delta0]; % Constant inputs for absolute model
U_p = @(t) U_k(t) - U_abs(t); % Perturbation input as deviation from nominal input
up = cell2mat(arrayfun(U_p, Tspan, 'UniformOutput', false));

%% Solve the Nonlinear Kinematic Bicycle Model using ODE45
[t_k, X_k] = ode45(@(t, X) kinematic_bicycle_model_2inputs(t, X, U_k(t), L), Tspan, X0_k);

%% Solve the Absolute State-Space Model using ODE45
[t_abs, X_abs] = ode45(@(t, X) absolute_ss_model_2inputs(t, X, U_abs(t), L), Tspan, X0_abs);

[t_p, X_p] = ode45(@(t, X) perturbation_ss_model_2inputs(t, X, U_p(t), L, ...
    interp1(t_abs, X_abs(:,4), t, 'linear', 'extrap'), ... % v_nom
    interp1(t_abs, X_abs(:,3), t, 'linear', 'extrap'), ... % psi_nom
    delta0), Tspan, X0_p);


%% Extract States
x_k = X_k(:,1); y_k = X_k(:,2); psi_k = X_k(:,3); v_k = X_k(:,4);
x_abs = X_abs(:,1); y_abs = X_abs(:,2); psi_abs = X_abs(:,3); v_abs = X_abs(:,4);
x_p = X_p(:,1); y_p = X_p(:,2); psi_p = X_p(:,3); v_p = X_p(:,4);

%% Plot Results
figure; hold on;
plot(x_k, y_k, 'b', 'LineWidth', 2);
plot(x_abs + x_p, y_abs + y_p, 'r--', 'LineWidth', 2);
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Kinematic', 'Nominal + Perturbation');
grid on;

%% Time Evolution of States
figure; hold on;
subplot(4,1,1); 
plot(t_k, v_k, 'k', 'LineWidth', 2); % 🔹 Corrected velocity plot
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Velocity Profile');
grid on;

subplot(4,1,2); 
plot(t_k, x_k, 'b', 'LineWidth', 2);
hold on; plot(t_p, x_p + x_abs, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)'); ylabel('X-position (m)');
legend('Kinematic', 'Nominal + Perturbation');

subplot(4,1,3); 
plot(t_k, y_k, 'b', 'LineWidth', 2);
hold on; plot(t_p, y_p + y_abs, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)'); ylabel('Y-position (m)');
legend('Kinematic', 'Nominal + Perturbation');

subplot(4,1,4); 
plot(t_k, psi_k, 'b', 'LineWidth', 2);
hold on; plot(t_p, psi_p + psi_abs, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)'); ylabel('Yaw Angle (Rad)');
legend('Kinematic', 'Nominal + Perturbation');
