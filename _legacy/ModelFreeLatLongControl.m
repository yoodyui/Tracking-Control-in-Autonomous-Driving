%This code demonstrates a model-free lateral control (PID/Stanley/Pure Pursuit)
% and model-free longitudinal control (PID) 
% implemented in the linearized model

clear; clc; close all;
% Parameters
L = 1.08;  % Wheelbase
N = 2000;
del_max = deg2rad(30);
Ts = 0.1;

% Initial state and input
u = [0.1*ones(1,N); 0*ones(1,N)];  % [a; delta (steering)]
x = [0; 0; 0; 0]; % [x; y; phi (heading);  v]

% Reference trajectory
n = 100; % Number of waypoints
[WP, velocities, headings] = waypoints_with_velocity_and_heading(n);
reference = [WP; zeros(1,length(WP))];

%Output matrix
C = [0 0 0 1];

%PID Gains
%Kp_lon = 0.9; Ki_lon = 0; Kd_lon = 0.1;
Kp = 0.1; Ki = 0; Kd = 0.02; Kt = 0.5; 
% Initialize Integral Term
integral_error = 0;
k_s = 0.1; %Stanley gain
%Pure Pursuit Parameters
Ld = 15; % Look-ahead distance in meters (larger values are suitable for high speeds)


for k = 1:N

    %True CTE
    %Find the closest waypoints
    [next(:,k),pass(:,k),next_index(k),pass_index(k)] = nextWP(WP,[x(1,k); x(2,k)]);
    

    %Calculate CTE
    e(k) = cte(next(:,k),pass(:,k),[x(1,k); x(2,k)]);

    %Calculate heading error by normalized angle difference
    theta_e(k) = mod(headings(next_index(k)) - x(3,k) + pi, 2*pi) - pi;


% Discrete-Time Kinematic Model
    x(:,k+1) = x(:,k) + Ts * [x(4,k) * cos(x(3,k));
                                       x(4,k) * sin(x(3,k));
                                       (x(4,k)/L) * tan(u(2,k));
                                       u(1,k)];
    % %Limit speed to be less than 3 m/s
    % if x(4,k+1) >= 6 %Change to 6 m/s to see the error growth
    %     x(4,k+1) = x(4,k);
    % end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1) Lateral Control
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Choice#1: PID controller (works best for all speeds)
    % integral_error = integral_error + e(k) * Ts; % Update Integral Term
    % if k==1
    %     u(2,k+1) = Kp*e(k) + Ki * integral_error + Kt*theta_e(k);
    % else
    %     u(2,k+1) = Kp*e(k) + Ki * integral_error + Kd * (e(k) - e(k-1)) / Ts + Kt * theta_e(k);
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Choice#2: Stanley Control Law (works better for low speeds less than 5m/s)
    %u(2,k+1) = theta_e(k) + atan2(k_s * e(k), u(1,k)); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Choice#3: Pure Pursuit Controller (works best for all speeds)
    % Compute the distances to all waypoints
    distances = sqrt((WP(1,:) - x(1,k)).^2 + (WP(2,:) - x(2,k)).^2);
    % Find the closest waypoint
    [~, closestIndex(k)] = min(distances);

    % Search for the look-ahead point starting from the closest waypoint
    searchIndex(k) = closestIndex(k);
    while true
        % Compute the distance from current position to the waypoint
        Ld_candidate(k) = sqrt((WP(1,searchIndex(k)) - x(1,k))^2 + (WP(2,searchIndex(k)) - x(2,k))^2);
        if Ld_candidate(k) >= Ld
            break;
        end
        searchIndex(k) = searchIndex(k) + 1;
        if searchIndex(k) > n
            searchIndex(k) = 1; % Wrap around
        end
        % Prevent infinite loop
        if searchIndex(k) == closestIndex(k)
            break;
        end
    end
    lookAheadIndex(k) = searchIndex(k);
    Ld_point(:,k) = WP(:, lookAheadIndex(k));

    %Compute Alpha (angle between vehicle heading and look-ahead point)
    delta_x = Ld_point(1,k) - x(1,k);
    delta_y = Ld_point(2,k) - x(2,k);
    angle_to_point = atan2(delta_y, delta_x);
    alpha(k) = angle_to_point - x(3,k);
    % Normalize alpha to [-pi, pi]
    alpha(k) = atan2(sin(alpha(k)), cos(alpha(k)));

    %Compute Steering Angle using Pure Pursuit Law
    u(2,k+1) = atan2(2 * L * sin(alpha(k)), Ld);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Limit -del_max < u(2,k+1) < del_max
    u(2,k+1) = max(min(u(2,k+1), del_max), -del_max);

     % Calculate distance to the initial position to terminate iteration
    distanceToStart(k) = dis([next(1,1); next(2,1)],[x(1,k); x(2,k)]);
    if distanceToStart(k) < 5 && k > 500
        fprintf('Stopping simulation at step %d because distance is less than the threshold.\n', k);
        break;
    end

    % % Plot next waypoint and current position
    % plot(next(1,k), next(2,k), 'bo'); hold on; % Next waypoint
    % plot(x(1,k), x(2,k), 'k*');               % Current position
    % 
    % % Add heading arrow
    % arrowLength = 1; % Length of the heading arrow
    % dx = arrowLength * cos(x(3,k)); % x-component of the heading direction
    % dy = arrowLength * sin(x(3,k)); % y-component of the heading direction
    % quiver(x(1,k), x(2,k), dx, dy, 'r', 'MaxHeadSize', 2); % Heading arrow
    % 
    % % Add heading arrow at each waypoint
    % dx = arrowLength * cos(headings(next_index(:,k))); % x-component of heading direction
    % dy = arrowLength * sin(headings(next_index(:,k))); % y-component of heading direction
    % quiver(WP(1, next_index(:,k)), WP(2, next_index(:,k)), dx, dy, 'g', 'MaxHeadSize', 2, 'AutoScale', 'off'); % Heading arrow
    % 
    % drawnow;

end


% Plot path tracking results
figure;
hold on;
plot(reference(1, :), reference(2, :), 'r--', 'LineWidth', 2); % Reference path
plot(x(1, :), x(2, :), 'b-.', 'LineWidth', 2); % Actual path followed
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Reference Path', 'Actual Path', 'Location', 'Best');
title('Path Tracking Performance');
grid on; box on; axis tight; 
hold off;

% Plot state variables over time
figure;

subplot(4, 1, 1);
plot(0:Ts:(length(x(1, :))-1)*Ts, x(1, :), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('X Position (m)');
title('X Position Over Time');
grid on; box on; axis tight; 

subplot(4, 1, 2);
plot(0:Ts:(length(x(2, :))-1)*Ts, x(2, :), 'g-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Y Position (m)');
title('Y Position Over Time');
grid on; box on; axis tight; 

subplot(4, 1, 3);
plot(0:Ts:(length(x(3, :))-1)*Ts, rad2deg(x(3, :)), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Heading (deg)');
title('Heading Angle Over Time');
grid on; box on; axis tight; 

subplot(4, 1, 4);
plot(0:Ts:(length(x(4, :))-1)*Ts, x(4, :), 'm-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Speed (m/s)');
title('Vehicle Speed Over Time');
grid on; box on; axis tight; 

% Plot control inputs (acceleration and steering angle)
figure;

subplot(2, 1, 1);
plot(0:Ts:(length(u(1, :))-1)*Ts, u(1, :), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Acceleration (m/s²)');
title('Acceleration Profile');
grid on; box on; axis tight; 

subplot(2, 1, 2);
plot(0:Ts:(length(u(2, :))-1)*Ts, rad2deg(u(2, :)), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Steering Angle (deg)');
title('Steering Angle Profile');
grid on; box on; axis tight; 

% Plot CTE and Heading Error Over Time
figure;

subplot(2, 1, 1);
plot(0:Ts:(length(e)-1)*Ts, e, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Cross-Track Error (m)');
title('Cross-Track Error Over Time');
grid on; box on; axis tight; 

subplot(2, 1, 2);
plot(0:Ts:(length(theta_e)-1)*Ts, rad2deg(theta_e), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Heading Error (deg)');
title('Heading Error Over Time');
grid on; box on; axis tight;  
