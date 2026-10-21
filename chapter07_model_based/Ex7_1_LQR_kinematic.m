% Example 7.1: Path-Indexed TV-LQR with a Kinematic Bicycle Model
%Lateral + Longitudinal control
clc; clear; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));
%Full control (Lateral+Longitudinal by LQR in Kinematic model)
%Successfully tracking if the control inputs are limit or not!!
%Picking the closest WP in the loop--Good for all velocities upto 100m/s
%nearest-waypoint / path-indexed TV-LQR

%% Vehicle Parameters
L = 2.8;           % Wheelbase (m)
%% Discrete-Time Simulation Parameters
Ts = 0.05;  % Sampling time (s)
T_final = 200;
N = floor(T_final / Ts);  % Number of steps
Tspan = (0:N-1) * Ts;  % Discrete time steps

%% Reference trajectory generation
[WP, velocities, headings] = waypoints_with_velocity_and_heading(N);
reference = [WP; unwrap(headings); velocities]; %unwrap reference heading

%% Initial Conditions
X_nom = reference; % Directly assign nominal states from the reference
X_lqr(:,1) = [WP(:,1); 0; 0]; % [x; y; yaw; velocity]
X_tilde = zeros(4,N);

%% Actuator limits
a_min = -3; a_max = 2;        % acceleration bounds [m/s^2]
delta_max = deg2rad(30);      % steering bound [rad]

%% Nominal inputs from the reference (path derivatives)
% The waypoints are spaced by distance, not by time, so the reference is
% differentiated along the arc length s:  d(.)/dt = v * d(.)/ds
s_ref   = [0, cumsum(vecnorm(diff(WP, 1, 2)))];          % arc length [m]
kappa   = gradient(reference(3,:), s_ref);               % curvature d(psi)/ds
dv_ds   = gradient(reference(4,:), s_ref);               % d(v)/ds
U_nom = zeros(2,N); %[acceleration;steering angle]
U_nom(1,:) = reference(4,:) .* dv_ds;                    % a_nom = v dv/ds
U_nom(2,:) = atan(L * kappa);                            % delta_nom = atan(L*kappa)
U_nom(1,:) = min(max(U_nom(1,:), a_min), a_max);         % keep the nominal feasible
U_nom(2,:) = min(max(U_nom(2,:), -delta_max), delta_max);

%% LQR Weight Matrices
%Q = diag([1, 1, 100, 10e-9]);
Q = diag([1, 1, 100, 10]);
R = diag([1, 1]);

%% Simulation Loop
for k = 1:N-1
    
    %% Compute lateral error from LQR
    [next_lqr(:,k),pass_lqr(:,k),next_index_lqr(k),pass_index_lqr(k)] = nextWP(WP,[X_lqr(1,k); X_lqr(2,k)]);
    e_lqr(k) = cte(next_lqr(:,k),pass_lqr(:,k),[X_lqr(1,k); X_lqr(2,k)]);
    %Calculate heading error from LQR
    theta_e_lqr(k) = mod(headings(next_index_lqr(k)) - X_lqr(3,k) + pi, 2*pi) - pi;

    %% Compute Longitudinal error from LQR along heading
    wp_error_vector = WP(:,next_index_lqr(k)) - [X_lqr(1,k); X_lqr(2,k)];
    long_error(k) = dot(wp_error_vector, [cos(X_lqr(3,k)); sin(X_lqr(3,k))]);

    % Update nominal states
    x = X_nom(1,next_index_lqr(k));
    y = X_nom(2,next_index_lqr(k));
    psi = X_nom(3,next_index_lqr(k));
    v_x = X_nom(4,next_index_lqr(k));

    %control inputs
    delta = U_nom(2,next_index_lqr(k));

    %% Linearization around nominal trajectory
    A = [1, 0, -Ts*v_x*sin(psi), Ts*cos(psi);
         0, 1,  Ts*v_x*cos(psi), Ts*sin(psi);
         0, 0,  1, Ts*tan(delta)/L;
         0, 0,  0, 1];

    B = [0, 0;
         0, 0;
         0, (Ts/L)*v_x*sec(delta)^2;
         Ts, 0];

    %% Controllability check (only initially)
    if k == 1 && rank(ctrb(A, B)) ~= size(A,1)
        error('The system is NOT controllable.');
    end

    %% LQR gain clearly computed
    K = dlqr(A, B, Q, R);

    %% Compute state deviation (error) from nominal
    X_tilde(1:2, k) = X_lqr(1:2,k) - X_nom(1:2,next_index_lqr(k));             % Position error
    X_tilde(3, k)   = wrapToPi(X_lqr(3,k) - X_nom(3,next_index_lqr(k)));       % Heading error 
    X_tilde(4, k)   = X_lqr(4,k) - X_nom(4,next_index_lqr(k));                 % Velocity error

    %% LQR input clearly computed
    U_lqr(:, k) = U_nom(:, next_index_lqr(k)) - K * X_tilde(:, k);

    %% Actuator limits (same bounds as in Chapter 6)
    U_lqr(1,k) = min(max(U_lqr(1,k), a_min), a_max);            % acceleration [m/s^2]
    U_lqr(2,k) = min(max(U_lqr(2,k), -delta_max), delta_max);   % steering [rad]

    %% Actual vehicle state update
    X_lqr(1,k+1) = X_lqr(1,k) + Ts * X_lqr(4,k) * cos(X_lqr(3,k));
    X_lqr(2,k+1) = X_lqr(2,k) + Ts * X_lqr(4,k) * sin(X_lqr(3,k));
    X_lqr(3,k+1) = X_lqr(3,k) + Ts * X_lqr(4,k)/L * tan(U_lqr(2,k));
    X_lqr(4,k+1) = X_lqr(4,k) + Ts * U_lqr(1,k);
    
    % Calculate distance to the initial position to terminate iteration
    distanceToStart(k) = dis([next_lqr(1,1); next_lqr(2,1)],[X_lqr(1,k); X_lqr(2,k)]);
    if distanceToStart(k) < 5 && k > 100
        fprintf('Stopping simulation at step %d because distance is less than the threshold.\n', k);
        break;
    end

    % % % Plot next waypoint and current position
    % plot(next_lqr(1,k), next_lqr(2,k), 'bo'); hold on; % Next waypoint
    % plot(X_lqr(1,k), X_lqr(2,k), 'k*');               % Current position
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
plot(X_lqr(1,:), X_lqr(2,:), 'r-.', 'LineWidth', 2);
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Reference Trajectory','LQR Trajectory');
grid on; box on;


%% Plot Control Inputs
figure; hold on;
subplot 211; hold on;
plot(Tspan(1:k), U_nom(1,next_index_lqr(1:k)));
plot(Tspan(1:k), U_lqr(1,1:k),'r--');
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
legend('Nominal','LQR'); grid on;box on; axis("tight");

subplot 212; hold on;
plot(Tspan(1:k), rad2deg(U_nom(2, next_index_lqr(1:k))));
plot(Tspan(1:k), rad2deg(U_lqr(2,1:k)),'r--');
xlabel('Time (s)'); ylabel('Steering Angle (deg)');
legend('Nominal','LQR'); grid on;box on; axis("tight");

%% Plot States
figure; hold on;
subplot 211; hold on;
plot(Tspan(1:k), rad2deg(X_nom(3,next_index_lqr(1:k))));
plot(Tspan(1:k), rad2deg(X_lqr(3,1:k)),'r--');
xlabel('Time (s)'); ylabel('Heading Angle (deg)');
legend('Nominal','LQR'); grid on; box on; axis("tight");

subplot 212; hold on;
plot(Tspan(1:k), X_nom(4,next_index_lqr(1:k)));
plot(Tspan(1:k), X_lqr(4,1:k),'r--');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Nominal','LQR'); grid on;box on; axis("tight");


%% Plot Errors
figure; hold on;
subplot 211;
plot(Tspan(1:k), e_lqr(1:k));
xlabel('Time (s)'); ylabel('CTE (m)');
grid on;  box on; axis("tight");
subplot 212; 
plot(Tspan(1:k), rad2deg(theta_e_lqr(1:k)));
xlabel('Time (s)'); ylabel('Heading Error (deg)');
grid on;  box on; axis("tight");

%% Save the four figures used in the book (figures/CH07)
figdir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'figures', 'CH07');
names  = {'FigEx6_1_path','FigEx6_1_control_input','FigEx6_1_states','FigEx6_1_errors'};
for i = 1:4
    exportgraphics(figure(i), fullfile(figdir, [names{i} '.jpeg']), 'Resolution', 300);
end
