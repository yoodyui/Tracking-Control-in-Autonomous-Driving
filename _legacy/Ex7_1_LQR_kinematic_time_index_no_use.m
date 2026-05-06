% Example 7.1 auxiliary script: time-indexed TV-LQR variant (not used in text)
clc; clear; close all;
%Full control (Lateral+Longitudinal by LQR in Kinematic model)
%Successfully tracking if the control inputs are limit or not!!
%Good for velocity below 10m/s
%time-indexed TV-LQR

%% Vehicle Parameters
L = 2.8%1.08;           % Wheelbase (m)
%% Discrete-Time Simulation Parameters
Ts = 0.01;  % Sampling time (s)
T_final = 200;
N = floor(T_final / Ts);  % Number of steps
Tspan = (0:N-1) * Ts;  % Discrete time steps

%% Reference trajectory generation
[WP, velocities, headings] = waypoints_with_velocity_and_heading(N);
reference = [WP; headings; velocities];

%% Initial Conditions
X_nom = reference; % Directly assign nominal states from the reference
X(:,1) = X_nom(:,1); % [x; y; yaw; velocity]
X_tilde = zeros(4,N);

%% Nominal inputs (compute from reference)
U_nom = zeros(2,N);
U_nom(1,1:N-1) = diff(reference(4,:))/Ts; % acceleration
U_nom(2,1:N-1) = diff(unwrap(reference(3,:)))./Ts .* (L ./ reference(4,1:end-1)); % steering angle
U_nom(1,N) = U_nom(1,N-1);
U_nom(2,N) = U_nom(2,N-1);

%% Find the nominal path
[t_k, X_k] = ode45(@(t, X) kinematic_bicycle_model_2inputs(t, X, ...
    [interp1(Tspan, U_nom(1,:), t, 'previous', 'extrap');...
    interp1(Tspan, U_nom(2,:), t, 'previous', 'extrap')],...
    L), Tspan, X(:,1));
X_k = X_k';

%% LQR Weight Matrices
%Q = diag([1, 1, 100, 10e-9]);
Q = diag([1, 1, 100, 10]);
R = diag([1, 1]);

%% Simulation Loop
for k = 1:N-1

    %% Linearization around nominal trajectory
    A = [1, 0, -Ts*X_nom(4,k)*sin(X_nom(3,k)), Ts*cos(X_nom(3,k));
         0, 1,  Ts*X_nom(4,k)*cos(X_nom(3,k)), Ts*sin(X_nom(3,k));
         0, 0,  1, Ts*tan(U_nom(2,k))/L;
         0, 0,  0, 1];

    B = [0, 0;
         0, 0;
         0, (Ts/L)*X_nom(4,k)*sec(U_nom(2,k))^2;
         Ts, 0];

    %% Controllability check (only initially)
    if k == 1 && rank(ctrb(A, B)) ~= size(A,1)
        error('The system is NOT controllable.');
    end

    %% LQR gain clearly computed
    K = dlqr(A, B, Q, R);

    %% Compute state deviation (error) from nominal
    %X_tilde(:, k) = X(:,k) - X_nom(:,k); %This wrapped angle creates little circular path at left top corner 
    X_tilde(1:2, k) = X(1:2,k) - X_nom(1:2,k);             % Position error
    X_tilde(3, k)   = wrapToPi(X(3,k) - X_nom(3,k));       % Heading error wrapped [-pi, pi]
    X_tilde(4, k)   = X(4,k) - X_nom(4,k);                 % Velocity error

    %% LQR input clearly computed
    U_lqr(:, k) = U_nom(:, k) - K * X_tilde(:, k);

    %% Apply saturation limits explicitly
    %U_lqr(1,k) = min(max(U_lqr(1,k), -3), 3); % acceleration limits
    %U_lqr(:,k) = min(max(U_lqr(:,k), deg2rad(-30)), deg2rad(30)); % steering limits

    %% Actual vehicle state update
    X(1,k+1) = X(1,k) + Ts * X(4,k) * cos(X(3,k));
    X(2,k+1) = X(2,k) + Ts * X(4,k) * sin(X(3,k));
    X(3,k+1) = wrapToPi(X(3,k) + Ts * X(4,k)/L * tan(U_lqr(2,k)));
    X(4,k+1) = X(4,k) + Ts * U_lqr(1,k);

    % % % Plot next waypoint and current position
    % plot(X_nom(1,k), X_nom(2,k), 'bo'); hold on; % Next waypoint
    % plot(X(1,k), X(2,k), 'k*');               % Current position
    % drawnow;


end

%% Plot Results: Trajectory Comparison
figure;
plot(reference(1,:), reference(2,:), 'b:', 'LineWidth', 2); hold on;
plot(X_k(1,:), X_k(2,:), 'r-.', 'LineWidth', 2);
plot(X(1,:), X(2,:), 'k:', 'LineWidth', 2);
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Reference Trajectory','Nominal Trajectory', 'LQR-Controlled Trajectory');
grid on; axis equal; box on;


%% Plot Control Inputs
figure; hold on;
subplot 211; hold on;
plot(Tspan(1:k), rad2deg(U_nom(1,1:k)));
plot(Tspan(1:k), rad2deg(U_lqr(1,1:k)),'r--');
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
legend('Nominal','LQR'); grid on;box on;

subplot 212; hold on;
plot(Tspan(1:k), rad2deg(U_nom(2,1:k)));
plot(Tspan(1:k), rad2deg(U_lqr(2,1:k)),'r--');
xlabel('Time (s)'); ylabel('Steering Angle (deg)');
legend('Nominal','LQR'); grid on;box on;

%% Plot States
figure; hold on;
subplot 211; hold on;
plot(Tspan(1:k), rad2deg(X_nom(3,1:k)));
plot(Tspan(1:k), rad2deg(X(3,1:k)),'r--');
xlabel('Time (s)'); ylabel('Heading (deg)');
legend('Nominal','LQR'); grid on; box on;

subplot 212; hold on;
plot(Tspan(1:k), X_nom(4,1:k));
plot(Tspan(1:k), X(4,1:k),'r--');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Nominal','LQR'); grid on;box on;
