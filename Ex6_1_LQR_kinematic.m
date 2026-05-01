%Lateral + Longitudinal control
clc; clear; close all;
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

%% Nominal inputs (compute from reference)
U_nom = zeros(2,N); %[acceleration;steering angle]
U_nom(1,1:N-1) = diff(reference(4,:))/Ts; % acceleration
U_nom(2,1:N-1) = diff(reference(3,:))./Ts .* (L ./ reference(4,1:end-1)); % steering angle
U_nom(1,N) = U_nom(1,N-1);
U_nom(2,N) = U_nom(2,N-1);

%% Find the nominal path
[t_k, X_k] = ode45(@(t, X) kinematic_bicycle_model_2inputs(t, X, ...
    [interp1(Tspan, U_nom(1,:), t, 'previous', 'extrap');...
    interp1(Tspan, U_nom(2,:), t, 'previous', 'extrap')],...
    L), Tspan, X_lqr(:,1));
X_k = X_k';

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

    %% Apply saturation limits explicitly
    %U_lqr(1,k) = min(max(U_lqr(1,k), -3), 3); % acceleration limits
    %U_lqr(:,k) = min(max(U_lqr(:,k), deg2rad(-30)), deg2rad(30)); % steering limits

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
%plot(X_k(1,:), X_k(2,:), 'm:', 'LineWidth', 2);
plot(X_lqr(1,:), X_lqr(2,:), 'r-.', 'LineWidth', 2);
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Reference Trajectory','LQR Trajectory');
grid on; box on;


%% Plot Control Inputs
figure; hold on;
subplot 211; hold on;
plot(Tspan(1:k), rad2deg(U_nom(1,next_index_lqr(1:k))));
plot(Tspan(1:k), rad2deg(U_lqr(1,1:k)),'r--');
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
