%Example3_3
%This code shows the path generated from the continuous time and discrete
%time Kinematic models are matched
clc; clear; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

%% Vehicle Parameters
L = 3.0;          % Wheelbase (m)
a_0 = 1.0;        % Nominal acceleration
psi_0 = 0.1;      % Nominal yaw angle (rad)
delta_0 = 5 * pi / 180;  % Nominal steering angle (rad)

%% Initial Conditions
X0_k = [0; 0; 0; 0]; % Initial state for kinematic model: [x; y; psi; v]
X0_tilde = [0; 0; 0; 0]; % Initial state for perturbation model

%% Discrete-Time Simulation Parameters
Ts = 0.1;         % Sampling time (s)
T_final = 50;     % Total simulation time (s)
N = floor(T_final / Ts);  % Number of discrete steps
Tspan = (0:N-1) * Ts;  % Discrete time steps

%% Define Input Signals
U_k = [a_0*sin(0.03*Tspan); delta_0 + (1*pi/180)*sin(0.5*Tspan)];  % Acceleration and steering angle inputs
U_nominal = [a_0 * ones(1, N); delta_0 * ones(1, N)]; % Nominal input for perturbation analysis

%% Initialize State Storage
X_k_d = zeros(4, N);
X_k_c = zeros(4, N);
X_tilde_d = zeros(4, N);
X_tilde_c = zeros(4, N);

%% Define State-Space Matrices (Nominal Model)
A = [0, 0, -a_0*sin(psi_0), cos(psi_0);
     0, 0,  a_0*cos(psi_0), sin(psi_0);
     0, 0,  0, (1/L) * tan(delta_0);
     0, 0,  0, 0];

B = [0, 0;
     0, 0;
     0, (a_0/L) * sec(delta_0)^2;
     1, 0];

%% Discretization (Euler Approximation)
A_d = eye(4) + Ts * A;
B_d = Ts * B;

%% Simulation Loop
for k = 1:N-1
    % Continuous-time Nonlinear Kinematic Model
    [~, X_k] = ode45(@(t, X) kinematic_bicycle_model_2inputs(t, X, U_k(:,k:k+1), L), Tspan(k:k+1), X_k_d(:,k));
    X_k_c(:,k+1) = X_k(end,:)';

    % Discrete-Time Kinematic Model
    X_k_d(:,k+1) = X_k_d(:,k) + Ts * [X_k_d(4,k) * cos(X_k_d(3,k));
                                       X_k_d(4,k) * sin(X_k_d(3,k));
                                       (X_k_d(4,k)/L) * tan(U_k(2,k));
                                       U_k(1,k)];

    % Continuous-Time Perturbation Model
    [~, X_tilde] = ode45(@(t, X) A * X + B * (U_k(:,k) - U_nominal(:,k)), Tspan(k:k+1), X_tilde_d(:,k));
    X_tilde_c(:,k+1) = X_tilde(end,:)';

    % Discrete-Time Perturbation Model
    X_tilde_d(:,k+1) = A_d * X_tilde_d(:,k) + B_d * (U_k(:,k) - U_nominal(:,k));
end

%% Plot Results
figure;
plot(X_k_d(1,:), X_k_d(2,:), 'b', 'LineWidth', 2); hold on;
plot(X_k_c(1,:), X_k_c(2,:), 'r--', 'LineWidth', 2);
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Discrete', 'Continuous');
grid on;

figure;
subplot(2,1,1);
plot(Tspan, X_k_d(3,:), 'b', 'LineWidth', 2); hold on;
plot(Tspan, X_k_c(3,:), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Yaw Angle (Rad)');
legend('Discrete', 'Continuous');
grid on;

subplot(2,1,2);
plot(Tspan, X_k_d(4,:), 'b', 'LineWidth', 2); hold on;
plot(Tspan, X_k_c(4,:), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Discrete', 'Continuous');
grid on;

%% Perturbation State Evolution
figure;
subplot(2,1,1);
plot(Tspan, X_tilde_d(1,:), 'b', 'LineWidth', 2); hold on;
plot(Tspan, X_tilde_c(1,:), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Perturbation in X Position');
legend('Discrete', 'Continuous');
grid on;

subplot(2,1,2);
plot(Tspan, X_tilde_d(3,:), 'b', 'LineWidth', 2); hold on;
plot(Tspan, X_tilde_c(3,:), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Perturbation in Yaw Angle');
legend('Discrete', 'Continuous');
grid on;

%% Error Analysis
figure;
subplot(2,1,1);
plot(Tspan, X_tilde_d(1,:) - X_tilde_c(1,:), 'k', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Error in X Position');
grid on;
title('Approximation Error in X Position');

subplot(2,1,2);
plot(Tspan, X_tilde_d(3,:) - X_tilde_c(3,:), 'k', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Error in Yaw Angle');
grid on;
title('Approximation Error in Yaw Angle');

disp('Simulation Complete');
