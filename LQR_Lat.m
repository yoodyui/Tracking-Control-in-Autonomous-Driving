%This code shows Lateral control by LQR with kinematic model--the nominal
%model has fixed input (steering angle).
%Not work!!!The system is uncontrollable.
clc; clear; close all;

%% Vehicle Parameters
L = 3.0;        % Wheelbase (m)
delta_0 = 10 * pi / 180;  % Initial steering angle

%% Initial Conditions
X0 = [0.1; 0.1; 1]; % [x, y, psi]
X_nom(:,1) = X0;
X(:,1) = X0;
X_tilde(:,1) = X0;

%% Discrete-Time Simulation Parameters
Ts = 0.05;  % Sampling time
N = floor(70 / Ts);  % Number of steps
Tspan = (0:N-1) * Ts;  % Time steps

%% Fixed Velocity for Lateral Control Only
v = 1; % Fixed velocity (constant)

%% Input Signals (Steering Control Only)
U_nom = repmat(delta_0, 1, N); % Steering reference
U = U_nom + 0.2 * sin(0.5*Tspan); % Disturbed steering input

%% Discrete-Time Simulation Loop
for k = 1:N-1

    %% LQR Design (Updated for 3-State System)
    A = [1, 0, -Ts * v * sin(X_nom(3,k));
         0, 1,  Ts * v * cos(X_nom(3,k));
         0, 0,  1];

    B = [0;
         0;
         (Ts/L) * sec(U_nom(k))^2 * v];

    % Check Controllability
    if rank(ctrb(A, B)) ~= size(A,1)
        error('The system is NOT controllable.');
    end

    K(k,:) = dlqr(A, B, diag([10, 10, 1]), 1); % Lateral-only LQR

    % Compute LQR control input (Steering Only)
    U_k(k) = U(k) - K(k,:) * X_tilde(:, k);

    % Apply limits to steering input
    U_k(k) = min(U_k(k), 30*pi/180);
    U_k(k) = max(U_k(k), -30*pi/180);

    % Update nominal model (Lateral motion only)
    X_nom(1,k+1) = X_nom(1,k) + Ts * v * cos(X_nom(3,k));
    X_nom(2,k+1) = X_nom(2,k) + Ts * v * sin(X_nom(3,k));
    X_nom(3,k+1) = X_nom(3,k) + Ts/L * v * tan(U_nom(k));

    % Update actual model (LQR-controlled lateral motion)
    X(1,k+1) = X(1,k) + Ts * v * cos(X(3,k));
    X(2,k+1) = X(2,k) + Ts * v * sin(X(3,k));
    X(3,k+1) = X(3,k) + Ts/L * v * tan(U_k(k));

    % Update perturbation state X_tilde
    X_tilde(:, k+1) = X(1:3, k+1) - X_nom(1:3, k+1);
end

%% Plot Results: Trajectory Comparison
figure; hold on;
plot(X_nom(1,:), X_nom(2,:), 'b', 'LineWidth', 2);
plot(X(1,:), X(2,:), 'r--', 'LineWidth', 2);
xlabel('X Position (m)', 'FontSize', 12);
ylabel('Y Position (m)', 'FontSize', 12);
title('Lateral Control: Nominal vs. LQR-Controlled', 'FontSize', 14);
legend('Nominal Trajectory', 'Actual Trajectory (LQR-Controlled)', 'Location', 'best');
grid on;

%% Plot Control Input: Steering Over Time
figure;
plot(U_k * (180/pi), 'r', 'LineWidth', 2); % Convert rad to degrees
xlabel('Time Step (k)', 'FontSize', 12);
ylabel('Steering Angle (degrees)', 'FontSize', 12);
title('LQR-Controlled Steering Input', 'FontSize', 14);
grid on;
legend('\delta_k (Steering Angle)', 'Location', 'best');
