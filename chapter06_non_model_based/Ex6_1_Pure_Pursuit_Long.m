% Example 6.1: Pure Pursuit with Longitudinal PID Speed Control
clear; clc; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

% Parameters
L = 1.08;              % Wheelbase
N = 1000;               % Simulation steps
del_max = deg2rad(30); % Max steering angle
Ts = 0.05;             % Sampling time

% Reference trajectory
n = 100; % Number of waypoints
[WP, velocities, headings] = waypoints_with_velocity_and_heading(n);
reference = [WP; unwrap(headings); velocities];

% Initial state and input
u = [zeros(1,N+1); zeros(1,N+1)];    % [acceleration; steering] (preallocate N+1)
x = zeros(4, N+1);                   % [x; y; psi; v]
x(:,1) = [WP(1:2,1); 0; 0];

% --- Longitudinal PID (speed) ---
KvP = 0.8; KvI = 0.2; KvD = 0.0;     % example gains; tune as needed
a_min = -3.0;                        % m/s^2 (braking)
a_max =  2.0;                        % m/s^2 (throttle)
iv = 0;                              % integrator state
ev_prev = 0;                         % prev speed error for derivative

% Pure Pursuit
Ld = 10;                              % Look-ahead distance [m]

% Storage (for plots)
e     = zeros(1,N);
theta_e = zeros(1,N);
ev    = zeros(1,N);                   % speed error log
v_ref_log = zeros(1,N);               % desired speed used

for k = 1:N
    % Nearest waypoint for geometric refs
    [next(:,k), pass(:,k), next_index(k), pass_index(k)] = nextWP(WP, [x(1,k); x(2,k)]);
    
    % Desired states from nearest index
    x_desired(1,k) = reference(1, next_index(k));  % X_ref
    x_desired(2,k) = reference(2, next_index(k));  % Y_ref
    x_desired(3,k) = reference(3, next_index(k));  % psi_ref
    x_desired(4,k) = reference(4, next_index(k));  % v_ref
    v_ref = x_desired(4,k);
    v_ref_log(k) = v_ref;

    % Compute Cross-Track Error (CTE) and Heading Error
    e(k) = cte(next(:,k), pass(:,k), [x(1,k); x(2,k)]);
    theta_e(k) = mod(x_desired(3,k) - x(3,k) + pi, 2*pi) - pi;

    % --- Longitudinal PID (speed control) ---
    ev(k) = v_ref - x(4,k);                   % speed error
    dev   = (ev(k) - ev_prev)/Ts;             % derivative
    iv_try = iv + ev(k)*Ts;                   % tentative integrator
    
    a_unsat = KvP*ev(k) + KvI*iv_try + KvD*dev;

    % Saturate and anti-windup (simple conditional integration)
    a_sat = min(max(a_unsat, a_min), a_max);
    if a_unsat == a_sat
        iv = iv_try;                  % accept integrator only if not saturating
    end
    u(1,k) = a_sat;                   % use current accel command
    ev_prev = ev(k);

    % --- Lateral: Pure Pursuit steering ---
    % Find look-ahead point ~ Ld ahead along path
    distances = hypot(WP(1,:) - x(1,k), WP(2,:) - x(2,k));
    [~, closestIndex] = min(distances);
    searchIndex = closestIndex;
    while true
        Ld_candidate = hypot(WP(1,searchIndex) - x(1,k), WP(2,searchIndex) - x(2,k));
        if Ld_candidate >= Ld, break; end
        searchIndex = searchIndex + 1;
        if searchIndex > n, searchIndex = 1; end
        if searchIndex == closestIndex, break; end
    end
    lookAheadIndex = searchIndex;
    Ld_point = WP(:, lookAheadIndex);

    % Alpha: angle between vehicle heading and target point
    dx = Ld_point(1) - x(1,k);
    dy = Ld_point(2) - x(2,k);
    alpha = atan2(dy, dx) - x(3,k);
    alpha = atan2(sin(alpha), cos(alpha));  % wrap

    % Pure Pursuit law (steering for NEXT step, ZOH on u(2,k) in update)
    u(2,k) = atan2(2 * L * sin(alpha), Ld);
    u(2,k) = min(max(u(2,k), -del_max), del_max);

    % Errors for plots
    e(k) = cte(next(:,k), pass(:,k), x(1:2,k));
    theta_e(k) = atan2(sin(x_desired(3,k) - x(3,k)), cos(x_desired(3,k) - x(3,k)));

    % --- Discrete kinematic update (uses u(1,k) and u(2,k)) ---
    x(:,k+1) = x(:,k) + Ts * [ x(4,k)*cos(x(3,k));
                               x(4,k)*sin(x(3,k));
                               (x(4,k)/L)*tan(u(2,k));
                               u(1,k) ];

    % Early stop if back near start (optional)
    d0 = hypot(WP(1,1)-x(1,k), WP(2,1)-x(2,k));
    if d0 < 5 && k > 100
        fprintf('Stopping at k=%d (near start)\n',k); break;
    end

    % % Live markers (optional)
    % plot(WP(1,lookAheadIndex), WP(2,lookAheadIndex), 'bo'); hold on;
    % plot(x(1,k), x(2,k), 'k*'); drawnow;
end

k = min(k, N);   % in case of early stop

% %% Plot styling
% set(groot,'defaultFigureColor','w', ...
%           'defaultAxesFontName','Helvetica', ...
%           'defaultAxesFontSize',16, ...
%           'defaultAxesFontWeight','bold', ...
%           'defaultAxesLineWidth',1.5, ...
%           'defaultLineLineWidth',2.5, ...
%           'defaultLegendFontSize',14, ...
%           'defaultLegendFontWeight','bold');

%% Trajectory
figure(1); hold on;
plot(reference(1,:), reference(2,:), 'b', 'LineWidth', 2); hold on;
plot(x(1,1:k), x(2,1:k), 'r--', 'LineWidth', 2);
%plot(x(1,1:k), x(2,1:k), '-.', 'LineWidth', 2);
xlabel('X (m)'); ylabel('Y (m)');
legend('Reference','Pure Pursuit + speed PID'); grid on; box on; axis equal;

%% States
t = (0:k-1)*Ts;
figure(2); hold on; 
subplot(2,1,1); hold on;
plot(t, rad2deg(x_desired(3,1:k)), 'b--');
plot(t, rad2deg(x(3,1:k)), 'r:');
%plot(t, rad2deg(x(3,1:k)), '-.');
xlabel('Time (s)'); ylabel('Heading (deg)');
legend('Reference','Actual'); grid on; box on; axis tight;

subplot(2,1,2); hold on;
plot(t, v_ref_log(1:k), 'b--');
plot(t, x(4,1:k), 'r:');
%plot(t, x(4,1:k), '-.');
xlabel('Time (s)'); ylabel('Speed (m/s)');
legend('Reference','Actual'); grid on; box on; axis tight;

%% Errors
figure(3); hold on;
subplot(3,1,1); hold on; plot(t, e(1:k)); ylabel('CTE (m)'); grid on; box on; axis tight;
subplot(3,1,2); hold on; plot(t, rad2deg(theta_e(1:k))); ylabel('Heading Err (deg)'); grid on; box on; axis tight;
subplot(3,1,3); hold on; plot(t, ev(1:k)); ylabel('Speed Err (m/s)'); xlabel('Time (s)'); grid on; box on; axis tight;

%% Inputs
figure(4); hold on;
subplot(2,1,1); plot(t, rad2deg(u(2,1:k))); ylabel('Steering Angle(deg)'); grid on; box on; axis tight;
subplot(2,1,2); plot(t, u(1,1:k)); ylabel('Acceleration (m/s^2)'); xlabel('Time (s)'); grid on; box on; axis tight;
