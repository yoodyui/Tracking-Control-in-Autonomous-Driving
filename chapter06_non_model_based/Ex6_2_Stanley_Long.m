% Example 6.2: Stanley Controller + PID Speed Loop (no feedforward)
clear; clc; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

%% Parameters
L       = 1.08;                 % Wheelbase
N       = 1000;                 % Simulation steps
del_max = deg2rad(30);          % Max steering angle
Ts      = 0.05;                 % Sampling time

%% Reference trajectory
n = 100;                        % Number of waypoints
[WP, velocities, headings] = waypoints_with_velocity_and_heading(n);
reference = [WP; unwrap(headings); velocities];

%% Initial state and inputs
x = [WP(1:2,1); 0; 0];          % [x; y; psi; v]
u = zeros(2, N+1);              % [accel; steer]
u(1,:) = 0;                     % accel starts at 0 (let PID command it)
u(2,:) = 0;

%% Stanley lateral parameters
k_s     = 0.5;                  % lateral gain
epsilon = 0.1;                  % small number to avoid division by zero

%% PID speed controller (tracks x_desired(4,k))
Kpv = 0.8; Kiv = 0.2; Kdv = 0.0;   % speed PID gains (tune as needed)
a_max = 2.0;                        % accel saturation [m/s^2]
a_min = -3.0;                       % brake (decel) saturation [m/s^2]
int_ev = 0;                         % speed error integral
prev_ev = 0;

%% Preallocate logs
x_desired = zeros(4, N);
e_cte   = zeros(1,N);
theta_e = zeros(1,N);

%% Simulation loop
for k = 1:N

    % Nearest waypoint (for references & CTE)
    [next(:,k), pass(:,k), next_index(k), pass_index(k)] = nextWP(WP, [x(1,k); x(2,k)]);

    % Desired states from nearest waypoint
    x_desired(1,k) = reference(1, next_index(k));    % x_ref
    x_desired(2,k) = reference(2, next_index(k));    % y_ref
    x_desired(3,k) = reference(3, next_index(k));    % psi_ref (unwrap'ed)
    x_desired(4,k) = reference(4, next_index(k));    % v_ref

    % Errors for lateral control
    e_cte(k)   = cte(next(:,k), pass(:,k), [x(1,k); x(2,k)]);
    theta_e(k) = wrapToPi(x_desired(3,k) - x(3,k));

    %------------------------%
    % 1) Speed PID
    %------------------------%
    ev(k)   = x_desired(4,k) - x(4,k);     % speed error
    dev  = (k==1) * 0 + (k>1) * (ev(k) - prev_ev)/Ts;

    % (simple anti-windup: only integrate when not saturating or when integral helps)
    a_cmd_unsat = Kpv*ev(k) + Kiv*(int_ev) + Kdv*dev;
    u(1,k) = min(max(a_cmd_unsat, a_min), a_max); % commanded acceleration

    % Anti-windup back-calculation 
    if abs(a_cmd_unsat - u(1,k)) < 1e-12
        int_ev = int_ev + ev(k)*Ts;   % integrate only if not saturated
    end
    prev_ev = ev(k);

    %------------------------%
    % 2) Stanley lateral law
    %------------------------%
    v_for_stanley = max(x(4,k), 0.01);  % avoid divide-by-zero
    u(2,k) = theta_e(k) + atan2(k_s * e_cte(k), v_for_stanley + epsilon);
    u(2,k) = max(min(u(2,k), del_max), -del_max);

    %------------------------%
    % 3) Kinematic update
    %------------------------%
    x(:,k+1) = x(:,k) + Ts * [ x(4,k)*cos(x(3,k));
                               x(4,k)*sin(x(3,k));
                               (x(4,k)/L)*tan(u(2,k));
                               u(1,k) ];

    % Early stop if we loop back near the start
    if k > 100
        d0 = hypot(x(1,k) - WP(1,1), x(2,k) - WP(2,1));
        if d0 < 5
            fprintf('Stopping at k=%d (near start)\n',k); 
            break;
        end
    end

    % Live markers (optional)
    % plot(x_desired(1,k), x_desired(2,k), 'bo'); hold on;
    % plot(x(1,k), x(2,k), 'k*'); drawnow;
end

k = min(k, N);  % in case we broke early

%% Plot styling
set(groot,'defaultFigureColor','w', ...
          'defaultAxesFontName','Helvetica', ...
          'defaultAxesFontSize',16, ...
          'defaultAxesFontWeight','bold', ...
          'defaultAxesLineWidth',1.5, ...
          'defaultLineLineWidth',3, ...
          'defaultLegendFontSize',14, ...
          'defaultLegendFontWeight','bold');

%% Trajectory
figure(1); hold on;
plot(reference(1,:), reference(2,:), 'b', 'LineWidth', 2); hold on;
plot(x(1,1:k), x(2,1:k), 'r--', 'LineWidth', 2);
%plot(x(1,1:k), x(2,1:k), ':', 'LineWidth', 2);
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Reference','Stanley + speed PID');
grid on; box on;

%% States (Heading & Velocity)
t = (0:Ts:(k-1)*Ts);
figure(2); hold on; 
subplot(2,1,1); hold on;
plot(t, rad2deg(x_desired(3,1:k)), 'b--', 'LineWidth', 2);
plot(t, rad2deg(x(3,1:k)), 'r:',  'LineWidth', 2);
%plot(t, rad2deg(x(3,1:k)), ':',  'LineWidth', 2);
xlabel('Time (s)'); ylabel('Heading (deg)');
legend('Reference','Actual'); grid on; box on; axis tight;

subplot(2,1,2); hold on;
plot(t, x_desired(4,1:k), 'b--', 'LineWidth', 2);
plot(t, x(4,1:k), 'r:',  'LineWidth', 2);
%plot(t, x(4,1:k), ':',  'LineWidth', 2);
xlabel('Time (s)'); ylabel('Speed (m/s)');
legend('Reference','Actual'); grid on; box on; axis tight;

%% Errors (CTE & Heading)
figure(3); hold on;
subplot(3,1,1); hold on; plot(t, e_cte(1:k)); ylabel('CTE (m)'); grid on; box on; axis tight;
subplot(3,1,2); hold on; plot(t, rad2deg(theta_e(1:k))); ylabel('Heading Err (deg)'); grid on; box on; axis tight;
subplot(3,1,3); hold on; plot(t, ev(1:k)); ylabel('Speed Err (m/s)'); xlabel('Time (s)'); grid on; box on; axis tight;

%% Control Inputs
figure(4); hold on; 
subplot(2,1,1);
plot(t, rad2deg(u(2,1:k)), 'LineWidth', 1.8);
xlabel('Time (s)'); ylabel('Steering Angle(deg)'); grid on; box on; axis tight;

subplot(2,1,2);
plot(t, u(1,1:k), 'LineWidth', 1.8);
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)'); grid on; box on; axis tight;
