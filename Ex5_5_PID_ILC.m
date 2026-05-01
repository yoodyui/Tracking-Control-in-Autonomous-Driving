% Example 5.5: PID + ILC Controller for Iterative Patheta_e Tracking
% Notation: u_j(k) = nominal PID steering, uL_j(k) = learned feedforward,
%           delta_j(k) = u_j(k) + uL_j(k), e_j(k) = cross-track error.
clear; clc; close all;

%% Parameters
L      = 1.08;              % Wheelbase [m]
Ts      = 0.05;               % Sampling time [s]
N       = 500;               % Steps per traversal
n_iter  = 20;                 % Number of iterations (traversals)
del_max = deg2rad(30);       % Steering saturation [rad]

%% Reference (same patheta_e each iteration)
n_wp = 100;
[WP, velocities, headings] = waypoints_with_velocity_and_heading(n_wp);
reference = [WP; unwrap(headings); velocities];

%% Nominal PID (u_j(k))
Kp = 0.1;  Ki = 0.0;  Kd = 0;  Kt = 1;

%% ILC update: uL_{j+1}(k) = uL_j(k) + L_ilc * e_j(k)
L_ilc = 0.01;                           % learning gain

%% Storage across iterations
x_all = cell(1, n_iter);    % states per iteration
e_all = cell(1, n_iter);    % CTE per iteration

%% Initialize learned feedforward for iteration 1
uL = zeros(1, N);           % u^L_1(k) = 0

for j = 1:n_iter

    x = zeros(4, N+1);
    x(:,1) = [WP(1:2,1); 0; 0];  % [x; y; psi (heading); vx]
    
    % Inputs: u(1,k)=accel, u(2,k)=delta (applied steering)
    u = [1 * ones(1, N); zeros(1, N)];  % [acceleration; steering angle]
    
    integral_error = 0;
    e  = zeros(1, N);
    theta_e = zeros(1, N);   % heading error
    
    for k = 1:N
        % Nearest waypoint and segment for CTE computation
        [next(:, k), pass(:, k), next_index(k), pass_index(k)] = nextWP(WP, [x(1,k); x(2,k)]);

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


        % Nominal PID steering u_j(k)
        if k == 1
            u_pid = Kp*e(k) + Ki*integral_error + Kt*theta_e(k);
        else
            de = (e(k) - e(k-1)) / Ts;
            integral_error = integral_error + e(k)*Ts;
            u_pid = Kp*e(k) + Ki*integral_error + Kd*de + Kt*theta_e(k);
        end
        
        % Applied steering: delta_j(k) = u_j(k) + u^L_j(k)
        delta = u_pid + uL(k);
        delta = max(min(delta, del_max), -del_max);
        u(2,k) = delta;
        
        % Vehicle kinematics (bicycle model)
        x(:,k+1) = x(:,k) + Ts * [ x(4,k)*cos(x(3,k));
            x(4,k)*sin(x(3,k));
            (x(4,k)/L)*tan(delta);
            u(1,k) ];

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
    
    %close all;
    % Save iteration data
    x_all{j} = x;
    e_all{j} = e;
    rmse_cte(j) = sqrt(mean(e.^2));   % RMSE for iteration j

    % ILC update for next iteration
    uL = uL + L_ilc * e;
end

%% Plots: path tracking across iterations
figure; hold on; grid on; box on;
plot(WP(1,:), WP(2,:), 'k--', 'LineWidth', 2);
C = lines(n_iter);
for j = 1:n_iter
    plot(x_all{j}(1,:), x_all{j}(2,:), 'Color', C(j,:), 'LineWidth', 1.5);
end
xlabel('X (m)'); ylabel('Y (m)');
title('PID + ILC Iterative Path Tracking');
lgd = [{'Reference'}, arrayfun(@(jj) sprintf('Iter %d', jj), 1:n_iter, 'UniformOutput', false)];
legend(lgd{:}, 'Location', 'best');

%% 3D plot of CTE across iterations (line plot)
figure; hold on; grid on; box on;
for j = 1:n_iter
    t = (0:numel(e_all{j})-1)*Ts;       % time axis
    plot3(t, (j-1)*ones(size(t)), e_all{j}, 'LineWidth', 1.5);
end
xlabel('Time t (s)');
ylabel('Iteration j');
zlabel('Cross-Track Error e_j(k) (m)');
title('CTE Trajectories Across Iterations');
view(135, 25);   % adjust angles for better view


%% RMSE
figure; 
plot(0:n_iter-1, rmse_cte, '-o', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on; box on;
xlabel('Iteration j');
ylabel('RMSE(e_j) (m)');
title('CTE RMSE per Iteration (PID+ILC)');
