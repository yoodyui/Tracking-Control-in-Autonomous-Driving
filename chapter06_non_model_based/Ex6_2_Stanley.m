% Example 6.2: Stanley Controller Implementation for Lateral Control
clear; clc; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

% Parameters
L = 1.08;  % Wheelbase
N = 500; % Simulation steps
del_max = deg2rad(30); % Max steering angle
Ts = 0.05;  % Sampling time

% Reference trajectory
n = 100; % Number of waypoints
[WP, velocities, headings] = waypoints_with_velocity_and_heading(n);
reference = [WP; unwrap(headings); velocities];

% Initial state and input
u = [1 * ones(1, N); 0 * ones(1, N)];  % [acceleration; steering angle]
x = [WP(1:2,1); 0; 0];  % [x; y; psi (heading); vx]

% Stanley Controller Parameters
k_s = 0.5;  % Gain for lateral error correction
epsilon = 0.1; % Small constant to avoid division by zero

% Simulation loop
for k = 1:N
    % Find the closest waypoints
    [next(:, k), pass(:, k), next_index(k), pass_index(k)] = nextWP(WP, [x(1, k); x(2, k)]);
    
    % Store desired states
    x_desired(1,k) = reference(1, next_index(k));  % Desired X position
    x_desired(2,k) = reference(2, next_index(k));  % Desired Y position
    x_desired(3,k) = reference(3, next_index(k)); % Desired heading
    x_desired(4,k) = reference(4, next_index(k)); % Desired speed

    %Set speed to the desired speed
    x(4,k) = x_desired(4,k);

    % Compute Cross-Track Error (CTE) and Heading Error
    e(k) = cte(next(:, k), pass(:, k), [x(1, k); x(2, k)]);
    theta_e(k) = mod(x_desired(3,k) - x(3, k) + pi, 2*pi) - pi;

    % Stanley Control Law
    u(2, k) = theta_e(k) + atan2(k_s * e(k), x(4, k) + epsilon);

    % Limit steering angle within range -del_max < u(2, k+1) < del_max
    u(2, k) = max(min(u(2, k), del_max), -del_max);


    % Update Kinematic Model
    x(:, k+1) = x(:, k) + Ts * [x(4, k) * cos(x(3, k));
                                 x(4, k) * sin(x(3, k));
                                 (x(4, k) / L) * tan(u(2, k));
                                 u(1, k)];

    % Calculate distance to the initial position to terminate iteration
    distanceToStart(k) = dis([WP(1,1); WP(2,1)],[x(1,k); x(2,k)]);
    if distanceToStart(k) < 5 && k > 100
        fprintf('Stopping simulation at step %d because distance is less than the threshold.\n', k);
        break;
    end

    % % % Plot next waypoint and current position
    % plot(x_desired(1,k), x_desired(2,k), 'bo'); hold on; % Next waypoint
    % plot(x(1,k), x(2,k), 'k*');               % Current position
    % drawnow;

end

%% Set Global defaults--affects all figures
set(groot,'defaultFigureColor','w', ...
          'defaultAxesFontName','Helvetica', ...
          'defaultAxesFontSize',16, ...
          'defaultAxesFontWeight','bold', ...
          'defaultAxesLineWidth',1.5, ...
          'defaultLineLineWidth',3, ...
          'defaultLegendFontSize',14, ...
          'defaultLegendFontWeight','bold');

%% Plot Results: Trajectory Comparison
figure;
plot(reference(1,:), reference(2,:), 'b', 'LineWidth', 2); hold on;
plot(x(1,1:k), x(2,1:k), 'r--', 'LineWidth', 2);
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Reference Trajectory','Stanley Trajectory');
grid on; box on; 

%% Plot States (Heading & Velocity)
figure; hold on;
subplot(2,1,1); hold on;
plot((0:Ts:(k-1)*Ts), rad2deg(x_desired(3,1:k)), 'b--', 'LineWidth', 2);
plot((0:Ts:(k-1)*Ts), rad2deg(x(3,1:k)), 'r:', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Heading Angle(deg)');
legend('Reference','Stanley','Location','best');
grid on; box on; axis("tight");

subplot(2,1,2); hold on;
plot((0:Ts:(k-1)*Ts), x_desired(4,1:k), 'b--', 'LineWidth', 2);
plot((0:Ts:(k-1)*Ts), x(4,1:k), 'r:', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Reference','Stanley','Location','best');
grid on; box on; axis("tight");


%% Plot Errors (CTE, Heading Error, Longitudinal Error)
figure; hold on;
subplot(2,1,1);
plot((0:Ts:(k-1)*Ts), e(1:k), 'LineWidth', 1.8);
xlabel('Time (s)'); ylabel('CTE (m)');
grid on; box on; axis("tight");

subplot(2,1,2);
plot((0:Ts:(k-1)*Ts), rad2deg(theta_e(1:k)), 'LineWidth', 1.8);
xlabel('Time (s)'); ylabel('Heading Error (deg)');
grid on; box on; axis("tight");


%% Plot Control Inputs (Steering)
figure; hold on;
plot((0:Ts:(k-1)*Ts), rad2deg(u(2,1:k)), 'LineWidth', 1.8);
xlabel('Time (s)'); ylabel('Steering Angle(deg)');
grid on; box on; axis("tight");
