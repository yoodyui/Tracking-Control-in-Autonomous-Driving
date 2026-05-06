% Example 6.4: PID Controller with Dynamic Bicycle Model
%Update the nearest WP to the current position in the loop
%Only Lat control (use this version should be enough)
%Good for velocity below 40m/s (>40m/s results in wobbling paths)
clc; clear; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

%% Vehicle Parameters
m = 1500;       % Vehicle mass (kg)
Lf = 1.2;       % Distance from CG to front axle (m)
Lr = 1.6;       % Distance from CG to rear axle (m)
L = Lf + Lr;
Iz = m * ((0.5 * (Lf + Lr))^2);  % Yaw moment of inertia (kg*m^2)
Cd = 0.32;      % Drag coefficient
A = 2.2;        % Frontal area (m^2)
rho = 1.225;    % Air density (kg/m^3)
CxA = 0.5 * Cd * A * rho;
Cxf = 60000;
Cxr = 120000;
Cyf = 30000;
Cyr = 30000;

%% Simulation Time
Ts = 0.05;
T_final = 200;
N = floor(T_final / Ts);
Tspan = (0:N-1) * Ts;

%% Reference trajectory generation
[WP, velocities, headings] = waypoints_with_velocity_and_heading(N);
reference = [WP; headings; velocities];

%% Initialize states
X = zeros(6, N);
%X(:,1) = [WP(:,1); headings(1); velocities(1); 0; 0];

%% PID Controller Parameters
Kp = 0.1; Ki = 0; Kd = 0.02; Kt = 1; 
integral_error = 0;

%% Simulation Loop
for k = 1:N-1
    % Current state
    x = X(1,k);
    y = X(2,k);
    psi = X(3,k);
    v_x = X(4,k);
    v_y = X(5,k);
    omega = X(6,k);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1) Lateral Control
    %% Compute lateral error
    [next(:,k),pass(:,k),next_index(k),pass_index(k)] = nextWP(WP,[x; y]);
    e(k) = cte(next(:,k),pass(:,k),[x; y]);
    %Calculate heading error by normalized angle difference
    theta_e(k) = mod(headings(next_index(k)) - psi + pi, 2*pi) - pi;

    %% Longitudinal Distance Error (along heading)
    wp_error_vector = WP(:,next_index(k)) - [x; y];
    long_error(k) = dot(wp_error_vector, [cos(psi); sin(psi)]);

    %PID controller
    integral_error = integral_error + e(k) * Ts; % Update Integral Term

    if k==1
        delta = Kp*e(k) + Ki * integral_error + Kt*theta_e(k);
    else
        delta = Kp*e(k) + Ki * integral_error + Kd * (e(k) - e(k-1)) / Ts + Kt * theta_e(k);
    end

    % PID controller steering
    delta = min(max(delta, deg2rad(-30)), deg2rad(30));

    %Compute Sr dynamically based on desired acceleration
    accel_command = (velocities(next_index(k)) - v_x) / Ts;

    % Compute rear slip ratio Sr (RWD vehicle)
    Sr(k) = (m * accel_command) / (2 * Cxr);

    % Front slip Sf = 0 for RWD
    Sf = 0;

    % Assume zero slip ratios for simplicity
    U(:,k) = [delta; Sf; Sr(k)];


    % Dynamic bicycle model update
    dX = dynamic_bicycle_model(0, X(:,k), U(:,k), m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr);
    X(:,k+1) = X(:,k) + Ts * dX;


    % Calculate distance to the initial position to terminate iteration
    distanceToStart(k) = dis([next(1,1); next(2,1)],[x; y]);
    if distanceToStart(k) < 5 && k > 100
        fprintf('Stopping simulation at step %d because distance is less than the threshold.\n', k);
        break;
    end

    % % Plot next waypoint and current position
    % plot(next(1,k), next(2,k), 'bo'); hold on; % Next waypoint
    % plot(X(1,k), X(2,k), 'k*');               % Current position
    % drawnow;

end

%% Plot Trajectory
figure;
plot(reference(1,:), reference(2,:), 'b:', 'LineWidth', 2); hold on;
plot(X(1,1:k), X(2,1:k), 'r-.', 'LineWidth', 2);
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Reference Trajectory', 'PID Controlled Trajectory');
grid on;
axis equal;

%% Plot States
figure;
subplot(2,1,1); hold on;
plot(Tspan(1:k), unwrap(rad2deg(headings(next_index(1:k))))); xlim('tight'); 
plot(Tspan(1:k), unwrap(rad2deg(X(3,1:k))),'r--');
xlabel('Time (s)'); ylabel('Heading Angle (deg)');
legend('Desired','Actual'); grid on; box on;

subplot(2,1,2); hold on;
plot(Tspan(1:k), velocities(next_index(1:k))); xlim('tight'); 
plot(Tspan(1:k), X(4,1:k),'r--');
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)');
legend('Desired','Actual'); grid on; box on;

%% Plot Errors
figure;
subplot(211); plot(Tspan(1:k), e(1:k)); xlim('tight');   
xlabel('Time (s)'); ylabel('CTE (m)'); grid on;
subplot(212); plot(Tspan(1:k), rad2deg(theta_e(1:k))); xlim('tight'); 
xlabel('Time (s)'); ylabel('Heading Error (deg)'); grid on;

%% Plot Control Inputs
% figure; hold on;
% subplot 311; plot(Tspan(1:k), rad2deg(U(1,1:k)));
% xlabel('Time (s)'); ylabel('Steering Angle (deg)'); grid on; box on;
% subplot 312; plot(Tspan(1:k), rad2deg(U(2,1:k)));
% xlabel('Time (s)'); ylabel('Sf'); grid on; box on;
% subplot 313; plot(Tspan(1:k), rad2deg(U(3,1:k)));
% xlabel('Time (s)'); ylabel('Sr'); grid on; box on;
