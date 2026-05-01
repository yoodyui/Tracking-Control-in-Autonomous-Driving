%For Figures in Example 1, chapter 3
%The nonlinear kinematic model  with 3 states and 1 input can be splitted into two models: absolute
%and perturbation where the input to absolute is delta_0 and that for
%perturbation is u_k-u_abs (kinematic = absolute + perturbation)
%This code needs the following function files: 
% (1) kinematic_bicycle_model.m
% (2) absolute_ss_model.m
% (3) perturbation_ss_model.m

clc; clear; close all;

%% Vehicle Parameters
L = 3.0;        % Wheelbase (m)
v_0 = 10.0;        % Constant longitudinal velocity (m/s)
psi_0 = 0.1;    % Nominal yaw angle (rad)
delta_0 = 5 * pi / 180;  % Nominal steering angle (rad)

%% Initial Conditions
X0_k = [0; 0; 0]; % Initial state for kinematic model: [x; y; psi]
X0_abs = [0; 0; 0];   % Initial state for absolute model: [x; y; psi]
X0_p = [0; 0; 0]; % Initial state for perturbation model: [x; y; psi]

%% Simulation Time
Tspan = 0:0.1:45;  % Start and end time

%% Input Signal (Steering Angle)
U_k = @(t) (delta_0 + 1*pi/180*sin(0.5*t)); % Steering input for kinematic model (combination of absolute and perturbation)
U_abs = @(t) (delta_0); % Steering input for absolute model
U_p = @(t) (U_k(t) - U_abs(t)); % Steering input for perturbation model

%% Velocity
v = @(t) (v_0*sin(0.03*t));
v_values = arrayfun(v, Tspan); 
%% Solve the Nonlinear Bicycle Model using ODE45
[t_k, X_k] = ode45(@(t, X) kinematic_bicycle_model(t, X, U_k(t), L, v(t)), Tspan, X0_k);
%Solve the Absolute State-Space Model using ODE45
[t_abs, X_abs] = ode45(@(t, X) absolute_ss_model(t, X, U_abs(t), L, v(t)), Tspan, X0_abs);
%Solve the Perturbation Model using ODE45
[t_p, X_p] = ode45(@(t, X) perturbation_ss_model(t, X, U_p(t), L, v(t),...
    interp1(t_abs, X_abs(:,3), t, 'linear', 'extrap'), ...
    delta_0), Tspan, X0_p);


%% Extract States
x_abs = X_abs(:,1);
y_abs = X_abs(:,2);
psi_abs = X_abs(:,3);

x_k = X_k(:,1);
y_k = X_k(:,2);
psi_k = X_k(:,3);

x_p = X_p(:,1);
y_p = X_p(:,2);
psi_p = X_p(:,3);


%% Plot Results
figure; hold on;
plot(x_k, y_k, 'b', 'LineWidth', 2);
plot(x_p+x_abs, y_p+y_abs, 'r--', 'LineWidth', 2);
xlabel('X Position (m)'); ylabel('Y Position (m)');
%title('Vehicle Trajectory');
legend('Kinematic', 'Nominal + Perturbation');
grid on;


figure; hold on;
subplot(4,1,1); 
plot(Tspan, v_values, 'k', 'LineWidth', 2); % 🔹 Plot velocity over time
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


% figure; hold on; 
% subplot 311; plot(t_k,x_k,'b', 'LineWidth', 2);
% hold on; plot(t_p, x_p+x_abs,'r--', 'LineWidth', 2); grid on;
% xlabel('Time (s)'); ylabel('X-position (m)');
% legend('Kinematic', 'Absolute + Perturbation');
% 
% subplot 312; plot(t_k,y_k,'b', 'LineWidth', 2);
% hold on; plot(t_p, y_p+y_abs,'r--', 'LineWidth', 2);  grid on;
% xlabel('Time (s)'); ylabel('Y-position (m)');
% legend('Kinematic', 'Absolute + Perturbation');
% 
% subplot 313; plot(t_k,psi_k,'b', 'LineWidth', 2);
% hold on; plot(t_p, psi_p+psi_abs, 'r--', 'LineWidth', 2); grid on;
% xlabel('Time (s)'); ylabel('Yaw Angle (Rad)');
% legend('Kinematic', 'Absolute + Perturbation');