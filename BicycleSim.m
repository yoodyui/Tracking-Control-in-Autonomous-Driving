%This code needs the following function files: 
% (1) kinematic_bicycle_model.m
% (2) dynamic_bicycle_model.m

clc; clear; close all;

%% Vehicle Parameters
m = 1500;       % Vehicle mass (kg)
Lf = 1.5;       % Distance from CoG to front axle (m)
Lr = 1.5;       % Distance from CoG to rear axle (m)
Iz = m * ((0.5 * (Lf + Lr))^2);      % Yaw moment of inertia (kg*m^2)
Cd = 0.32;      % Drag coefficient
A = 2.2;        % Frontal area (m^2)
rho = 1.225;    % Air density (kg/m^3)
CA = 0.5 * Cd * A * rho;  % Aerodynamic resistance coefficient
Cxf = 100000;    % Longitudinal stiffness front (N)
Cxr = 100000;    % Longitudinal stiffness rear (N)
Cyf = 30000;    % Lateral stiffness front (N/rad)
Cyr = 30000;    % Lateral stiffness rear (N/rad)

%% Initial Conditions
x0 = 0;        % Initial x position (m)
y0 = 0;        % Initial y position (m)
psi0 = 0;      % Initial yaw angle (rad)
vx0 = 5;      % Initial longitudinal velocity (m/s)
vy0 = 0;       % Initial lateral velocity (m/s)
omega0 = 0;    % Initial yaw rate (rad/s)

X0 = [x0; y0; psi0; vx0; vy0; omega0];  % Initial state vector

%% Simulation Time
T_final = 20;  % Simulation time (s)
Tspan = [0 T_final];

%% Control Inputs (Constant Inputs for Demonstration)
delta = 5 * pi / 180;  % Steering angle (rad)
Sf = 0;%0.005;             % Front slip ratio
Sr = 0;%0.005;             % Rear slip ratio
U = [delta; Sf; Sr];   % Control input vector

%% Solve Nonlinear ODEs using ODE45
[t, X] = ode45(@(t, X) dynamic_bicycle_model(t, X, U, m, Iz, Lf, Lr, CA, Cxf, Cxr, Cyf, Cyr), Tspan, X0);

%% Extract Results
x = X(:,1);  % X position
y = X(:,2);  % Y position
psi = X(:,3); % Yaw angle
vx = X(:,4);  % Longitudinal velocity
vy = X(:,5);  % Lateral velocity
omega = X(:,6); % Yaw rate


%% Solve Kinematic Model using ODE45
[t_k, X_k] = ode45(@(t, X) kinematic_bicycle_model(t, X, delta, Lf+Lr, vx0), Tspan, [x0; y0; psi0]);

%% Extract Results
x_k = X_k(:,1);   % X position
y_k = X_k(:,2);   % Y position
psi_k = X_k(:,3); % Yaw angle
%% Plot Results
figure;
subplot(1,2,1);
plot(x_k, y_k, 'b--', 'LineWidth', 2); hold on;
plot(x, y, 'r:', 'LineWidth', 2);
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Vehicle Trajectory');
legend('Dynamic','Kinematic');
grid on;

subplot(1,2,2);
plot(t_k, psi_k * (180/pi), 'b--', 'LineWidth', 2); hold on;
plot(t, psi * (180/pi), 'r:', 'LineWidth', 2); 
xlabel('Time (s)'); ylabel('Yaw Angle (deg)');
title('Yaw Angle vs Time');
legend('Dynamic','Kinematic');
grid on;

figure;
subplot(3,1,1);
plot(t, vx, 'b--', 'LineWidth', 2); hold on;
plot(t, vx0*ones(size(vx)), 'r:', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)');
title('Longitudinal Velocity vs Time');
legend('Dynamic','Kinematic');
grid on;

subplot(3,1,2);
plot(t, vy, 'b--', 'LineWidth', 2); hold on;
plot(t, 0*vy, 'r:', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Lateral Velocity (m/s)');
title('Lateral Velocity vs Time');
legend('Dynamic','Kinematic');
grid on;

subplot(3,1,3);
plot(t, omega, 'b--', 'LineWidth', 2); hold on;
plot(t_k, [X_k(1,3); diff(X_k(:,3))./diff(t_k)], 'r:', 'LineWidth', 2)
xlabel('Time (s)'); ylabel('Yaw Rate (rad/s)');
title('Yaw Rate vs Time');
legend('Dynamic','Kinematic');
grid on;
