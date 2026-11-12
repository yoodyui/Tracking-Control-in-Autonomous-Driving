% Example 3.3: independent continuous-time and Euler trajectories.
% The linear perturbation model is frozen at the initial nominal point,
% as in the revised manuscript; it is not an LTV trajectory linearization.
clc; clear; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

%% Vehicle Parameters
L = 3.0;
v_0 = 10.0;
a_0 = 0.0;
psi_0 = 0.1;
delta_0 = 5 * pi / 180;

%% Initial Conditions
X0_k = [0; 0; psi_0; v_0];
X0_tilde = [0; 0; 0; 0];

%% Discrete-Time Simulation Parameters
Ts = 0.1;
T_final = 50;
N = floor(T_final / Ts);
Tspan = (0:N-1) * Ts;

%% Define Input Signals
U_k = [0.2*sin(0.03*Tspan); delta_0 + (1*pi/180)*sin(0.5*Tspan)];
U_nominal = [a_0 * ones(1, N); delta_0 * ones(1, N)];

%% Initialize State Storage
X_k_d = zeros(4, N);
X_k_c = zeros(4, N);
X_tilde_d = zeros(4, N);
X_tilde_c = zeros(4, N);
X_k_d(:,1) = X0_k;
X_k_c(:,1) = X0_k;
X_tilde_d(:,1) = X0_tilde;
X_tilde_c(:,1) = X0_tilde;

%% Define State-Space Matrices (Perturbation Model, frozen at the initial nominal point)
A = [0, 0, -v_0*sin(psi_0), cos(psi_0);
     0, 0,  v_0*cos(psi_0), sin(psi_0);
     0, 0,  0, (1/L) * tan(delta_0);
     0, 0,  0, 0];

B = [0, 0;
     0, 0;
     0, (v_0/L) * sec(delta_0)^2;
     1, 0];

%% Discretization (Euler Approximation)
A_d = eye(4) + Ts * A;
B_d = Ts * B;

%% Simulation Loop
for k = 1:N-1
    u = U_k(:,k);
    [~, X_k] = ode45(@(t, X) kinematic_bicycle_model_2inputs(t, X, u, L), ...
                     Tspan(k:k+1), X_k_c(:,k));
    X_k_c(:,k+1) = X_k(end,:)';

    % Forward-Euler nonlinear kinematic update
    x = X_k_d(:,k);
    X_k_d(:,k+1) = x + Ts * [x(4)*cos(x(3));
                             x(4)*sin(x(3));
                             (x(4)/L)*tan(u(2));
                             u(1)];

    u_tilde = u - U_nominal(:,k);
    [~, X_tilde] = ode45(@(t, X) A * X + B * u_tilde, ...
                         Tspan(k:k+1), X_tilde_c(:,k));
    X_tilde_c(:,k+1) = X_tilde(end,:)';
    X_tilde_d(:,k+1) = A_d * X_tilde_d(:,k) + B_d * u_tilde;
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
