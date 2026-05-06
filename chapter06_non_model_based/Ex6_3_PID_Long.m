% Example 6.3: PID Controller Implementation for Lateral + Longitudinal Control
clear; clc; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

% Parameters
L = 1.08;                 % Wheelbase
N = 1000;                  % Simulation steps
del_max = deg2rad(30);    % Max steering angle
Ts = 0.05;                % Sampling time

% Reference trajectory
n = 100; % Number of waypoints
[WP, velocities, headings] = waypoints_with_velocity_and_heading(n);
reference = [WP; unwrap(headings); velocities];

% Initial state and input
u = [zeros(1, N); zeros(1, N)];         % [acceleration; steering angle]
x = [WP(1:2,1); 0; 0];                  % [x; y; psi (heading); vx]

% --- Lateral PID (on cross-track) + heading P ---
Kp = 0.1;       % CTE P
Ki = 0.0;       % CTE I
Kd = 0.0;       % CTE D
Kt = 1.0;       % heading P
int_cte = 0;    % integral state for lateral error

% --- Longitudinal PID (speed control) ---
KvP = 0.8;      % speed P (tune)
KvI = 0.2;      % speed I (tune)
KvD = 0.05;     % speed D (tune)
a_min = -3.0;   % accel lower limit [m/s^2] (tune to platform)
a_max =  2.0;   % accel upper limit [m/s^2]
iv = 0;         % speed integrator
ev_prev = 0;    % previous speed error

% Simulation loop
for k = 1:N
    % Nearest waypoint (for geometric references)
    [next(:,k), pass(:,k), next_index(k), pass_index(k)] = nextWP(WP, [x(1,k); x(2,k)]);
    
    % Desired states (from nearest index)
    x_desired(1,k) = reference(1, next_index(k));   % X_ref
    x_desired(2,k) = reference(2, next_index(k));   % Y_ref
    x_desired(3,k) = reference(3, next_index(k));   % psi_ref
    x_desired(4,k) = reference(4, next_index(k));   % v_ref

    % --- Lateral errors ---
    e_cte(k) = cte(next(:,k), pass(:,k), [x(1,k); x(2,k)]);                       % CTE
    theta_e(k) = atan2( sin(x_desired(3,k) - x(3,k)), cos(x_desired(3,k) - x(3,k)) ); % heading error

    % --- Longitudinal PID (speed control) -> u(1,k) = acceleration ---
    v_ref = x_desired(4,k);
    ev(k)  = v_ref - x(4,k);       % speed error
    dev    = (k==1) * 0 + (k>1) * (ev(k) - ev_prev)/Ts;   % derivative
    iv_try = iv + ev(k)*Ts;                                 % tentative integrator
    a_unsat = KvP*ev(k) + KvI*iv_try + KvD*dev;
    % Saturation + simple anti-windup (clamp)
    a_sat = min(max(a_unsat, a_min), a_max);
    if abs(a_unsat - a_sat) < 1e-9
        iv = iv_try;      % accept integrator when not saturating
    end
    u(1,k) = a_sat;
    ev_prev = ev(k);

    % --- Lateral PID (steering) -> u(2,k) ---
    if k == 1
        u(2,k) = Kp*e_cte(k) + Ki*int_cte + Kt*theta_e(k);
    else
        de_cte = (e_cte(k) - e_cte(k-1))/Ts;
        int_cte = int_cte + e_cte(k)*Ts;
        u(2,k) = Kp*e_cte(k) + Ki*int_cte + Kd*de_cte + Kt*theta_e(k);
    end
    % Steering saturation
    u(2,k) = max(min(u(2,k), del_max), -del_max);

    % --- Discrete-Time Kinematic Model Update (uses u(:,k)) ---
    x(:, k+1) = x(:, k) + Ts * [ x(4, k) * cos(x(3, k));
                                 x(4, k) * sin(x(3, k));
                                 (x(4, k) / L) * tan(u(2, k));
                                 u(1, k) ];

    % Early exit if back near start (optional)
    distanceToStart(k) = dis([WP(1,1); WP(2,1)], [x(1,k); x(2,k)]);
    if distanceToStart(k) < 5 && k > 100
        fprintf('Stopping simulation at step %d because distance is less than the threshold.\n', k);
        break;
    end

    % Live markers (optional)
    % plot(x_desired(1,k), x_desired(2,k), 'bo'); hold on;
    % plot(x(1,k), x(2,k), 'k*'); drawnow;
end

% ========== Plotting ==========
set(groot,'defaultFigureColor','w', ...
          'defaultAxesFontName','Helvetica', ...
          'defaultAxesFontSize',16, ...
          'defaultAxesFontWeight','bold', ...
          'defaultAxesLineWidth',1.5, ...
          'defaultLineLineWidth',3, ...
          'defaultLegendFontSize',14, ...
          'defaultLegendFontWeight','bold');

% Trajectory
figure(1);
plot(reference(1,:), reference(2,:), 'b', 'LineWidth', 2); hold on;
plot(x(1,1:k), x(2,1:k), 'r--', 'LineWidth', 2);
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Reference Trajectory','PID Trajectory');
grid on; box on;

% States (heading & velocity)
t = (0:Ts:(k-1)*Ts);
figure(2); 
subplot(2,1,1); hold on;
plot(t, rad2deg(x_desired(3,1:k)), 'b', 'LineWidth', 2);
plot(t, rad2deg(x(3,1:k)), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Heading (deg)');
legend('Reference','PID','Location','best'); grid on; box on; axis tight;

subplot(2,1,2); hold on;
plot(t, x_desired(4,1:k), 'b', 'LineWidth', 2);
plot(t, x(4,1:k), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Reference','PID','Location','best'); grid on; box on; axis tight;

% Errors (CTE & heading)
figure(3);
subplot(3,1,1); plot(t, e_cte(1:k)); ylabel('CTE (m)'); grid on; box on; axis tight;
subplot(3,1,2); plot(t, rad2deg(theta_e(1:k))); ylabel('Heading Err (deg)'); grid on; box on; axis tight;
subplot(3,1,3); plot(t, ev(1:k)); ylabel('Speed Err (m/s)'); xlabel('Time (s)'); grid on; box on; axis tight;


% Inputs
figure(4); 
subplot(2,1,1); hold on;
plot(t, rad2deg(u(2,1:k)), 'LineWidth', 1.8);
xlabel('Time (s)'); ylabel('Steering (deg)'); grid on; box on; axis tight;

subplot(2,1,2); hold on;
plot(t, u(1,1:k), 'LineWidth', 1.8);
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)'); grid on; box on; axis tight;
