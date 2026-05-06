%This code shows Lateral and Longitudinal control by LQR--the nominal
%model has fixed inputs (acceleration and steering angle.

clc; clear; close all;

%% Vehicle Parameters
L = 3.0;        % Wheelbase (m)
a_0 = 0.05;
delta_0 = 5 * pi / 180;

%% Initial Conditions
X0 = [0; 0; 0; 1]; % Initial state for kinematic model: [x; y; psi; v]

%% Discrete-Time Simulation Parameters
Ts = 0.05;  % Sampling time (s)
N = floor(70 / Ts);  % Number of steps
Tspan = (0:N-1) * Ts;  % Discrete time steps

%% Initialize State Storage
X_nom(:,1) = X0;
X(:,1) = X0;
X_lqr(:,1) = X0;
X_tilde(:,1) = zeros(4, 1);

%% Input Signals
U_nom = repmat([a_0; delta_0], 1, N); % Reference control inputs
U(1,:) = U_nom(1,1) + 0.1*sin(0.3*Tspan); % Disturbed Acceleration
U(2,:) = U_nom(2,1) + 0.2*sin(0.5*Tspan); % Disturbed Steering Angle


%% Discrete-Time Simulation Loop
for k = 1:N-1


    %% LQR design
    A = [1, 0, -Ts * X_nom(4,k) * sin(X_nom(3,k)), Ts * cos(X_nom(3,k));
        0, 1,  Ts * X_nom(4,k) * cos(X_nom(3,k)), Ts * sin(X_nom(3,k));
        0, 0,  1, (Ts/L) * tan(U_nom(2,k));
        0, 0,  0, 1];

    B = [0, 0;
        0, 0;
        0, (Ts/L) * sec(U_nom(2,k))^2 * X_nom(4,k);
        Ts, 0];


    % Check Controllability
    if rank(ctrb(A, B)) ~= size(A,1)
        error('The system is NOT controllable.');
    end

    K(:,:,k) = dlqr(A, B, diag([10, 10, 1, 1]), diag([1, 1]));

    % Compute deviation in input
    %U_tilde(:,k) = U(:,k) - U_nom(:,k);

    % Compute state deviation using LQR
    U_k(:, k) = U(:, k) - K(:,:,k) * X_tilde(:, k);

    % Apply limits to control input
    U_k(:, k) = min(U_k(:, k), [1; 30*pi/180]);
    U_k(:, k) = max(U_k(:, k), [0; -30*pi/180]);

    % Update nominal model
    X_nom(1,k+1) = X_nom(1,k) + Ts * X_nom(4,k) * cos(X_nom(3,k));
    X_nom(2,k+1) = X_nom(2,k) + Ts * X_nom(4,k) * sin(X_nom(3,k));
    X_nom(3,k+1) = X_nom(3,k) + Ts/L * X_nom(4,k) * tan(U_nom(2,k));
    X_nom(4,k+1) = X_nom(4,k) + Ts * U_nom(1,k);

    % Update actual model
    X(1,k+1) = X(1,k) + Ts * X(4,k) * cos(X(3,k));
    X(2,k+1) = X(2,k) + Ts * X(4,k) * sin(X(3,k));
    X(3,k+1) = X(3,k) + Ts/L * X(4,k) * tan(U_k(2,k));
    X(4,k+1) = X(4,k) + Ts * U_k(1,k);

    % Update perturbation state X_tilde
    % X_tilde(:, k+1) = A * X_tilde(:, k) + B * U_tilde(:, k); %cannot use
    % this since it doesn't reflect full nonlinear system evolution
    X_tilde(:, k+1) = X(:, k+1) - X_nom(:, k+1);

end

%% Plot Results: Trajectory Comparison
figure; hold on;
plot(X_nom(1,:), X_nom(2,:), 'b', 'LineWidth', 2);
plot(X(1,:), X(2,:), 'r--', 'LineWidth', 2);
xlabel('X Position (m)', 'FontSize', 12);
ylabel('Y Position (m)', 'FontSize', 12);
title('Trajectory Comparison: Nominal vs. Actual', 'FontSize', 14);
legend('Nominal Trajectory', 'Actual Trajectory (LQR-Controlled)', 'Location', 'best');
grid on;

%Plot Results: Control Inputs Over Time
figure; 

% Plot Acceleration Input
subplot(2,1,1); hold on;
plot(U_k(1,:), 'b', 'LineWidth', 2);
xlabel('Time Step (k)', 'FontSize', 12);
ylabel('Acceleration Input (m/s^2)', 'FontSize', 12);
grid on;

% Plot Steering Input
subplot(2,1,2); hold on;
plot(U_k(2,:) * (180/pi), 'r', 'LineWidth', 2); % Convert rad to degrees
xlabel('Time Step (k)', 'FontSize', 12);
ylabel('Steering Angle (degrees)', 'FontSize', 12);
grid on;

