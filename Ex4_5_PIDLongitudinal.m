% Example: Discrete-Time PID Speed Tracking on a Straight Path (delta_k = 0)
% This script implements the control law described in:
%   a[k] = K_v^P e_v[k] + K_v^I sum_{i=0}^k e_v[i] Ts + K_v^D (e_v[k]-e_v[k-1])/Ts + K_ff * d(v_ref)/dt
% with acceleration saturation and conditional anti-windup.
%
% State (kinematic bicycle with zero steering):
%   x_{k+1} = x_k + Ts * v_k
%   y_{k+1} = y_k
%   phi_{k+1} = phi_k (== 0)
%   v_{k+1} = v_k + Ts * a_k_cmd

clear; clc; close all;

%% Parameters
L  = 1.08;                % Wheelbase [m] (unused when delta = 0)
N  = 300;                % Simulation steps
Ts = 0.1;                 % Sampling time [s]

% Acceleration limits
a_min = -0.5;%-3.0;             % [m/s^2] braking (lower bound)
a_max = 0.5;%2.0;             % [m/s^2] traction (upper bound)

% PID gains (exactly as notation in the chapter)
KvP = 0.9;                % K_v^P
KvI = 0.25;               % K_v^I
KvD = 0.05;               % K_v^D
Kff = 1.0;                % K_ff (set 0 to disable reference-rate feedforward)

%% Reference speed: trapezoidal profile (accelerate -> cruise -> decelerate)
t = (0:N-1)*Ts;
v_cruise = 6.0;           % [m/s] (~21.6 km/h)
alpha    = 0.6;           % [m/s^2] reference slope
t_rise   = v_cruise/alpha;
t_cruise = 10;            % [s] hold time at cruise speed
k1 = round(t_rise/Ts);
k2 = k1 + round(t_cruise/Ts);

vref = zeros(1,N);
for k = 1:N
    if k <= k1
        vref(k) = min(alpha*t(k), v_cruise);
    elseif k <= k2
        vref(k) = v_cruise;
    else
        vref(k) = max(v_cruise - alpha*(t(k)-t(k2)), 0);
    end
end

% Reference-rate (finite difference)
dvref = [0, diff(vref)/Ts];

%% Allocation / logs
% States: [x; y; phi; v], starting at rest, heading 0 rad, on x-axis
x = [0; 0; 0; 0];
x_log = zeros(4,N); x_log(:,1) = x;

% Control signals and errors
ev    = zeros(1,N);      % e_v[k] = v_ref[k] - v[k]
a_cmd = zeros(1,N);      % a[k] (after saturation)

% Integral accumulator for e_v (for K_v^I term)
Iev = 0;

%% Main simulation loop
for k = 1:N
    % Speed error at step k
    ev(k) = vref(k) - x(4);

    % Error derivative (finite difference)
    if k == 1
        dev = 0;
    else
        dev = (ev(k) - ev(k-1))/Ts;
    end

    % ------------------ PID law (exactly as in the chapter) ------------------
    % Pre-check unsaturated command with current integral term
    a_unsat_pre = KvP*ev(k) + KvI*Iev + KvD*dev + Kff*dvref(k);

    % % Conditional anti-windup gate:
    % % If command would saturate high and error is positive -> hold integral
    % % If command would saturate low  and error is negative -> hold integral
    % hold_int = (a_unsat_pre > a_max && ev(k) > 0) || ...
    %            (a_unsat_pre < a_min && ev(k) < 0);
    % 
    % if ~hold_int
        Iev = Iev + ev(k)*Ts;  % update integral sum
    % end

    % Final unsaturated command with possibly updated integral
    a_unsat = KvP*ev(k) + KvI*Iev + KvD*dev + Kff*dvref(k);

    % Saturation
    a_k = min(max(a_unsat, a_min), a_max);
    a_cmd(k) = a_k;
    % ------------------------------------------------------------------------

    % ------------------ Kinematic update (delta_k = 0) ----------------------
    % delta = 0 -> phi remains 0; y remains ~0
    x = x + Ts * [ x(4)*cos(x(3));    % x position
                   x(4)*sin(x(3));    % y position (stays 0)
                   (x(4)/L)*tan(0);   % heading (stays 0)
                   a_k ];             % speed update
    % ------------------------------------------------------------------------

    x_log(:,k) = x;
end

%% Plots
set(groot,'defaultFigureColor','w', ...
          'defaultAxesFontName','Helvetica', ...
          'defaultAxesFontSize',16, ...
          'defaultAxesFontWeight','bold', ...
          'defaultAxesLineWidth',1.5, ...
          'defaultLineLineWidth',3, ...
          'defaultLegendFontSize',14, ...
          'defaultLegendFontWeight','bold');

% Speed tracking
figure('Name','Speed Tracking','Color','w');
plot(t, vref, 'k--', 'LineWidth', 1.4); hold on;
plot(t, x_log(4,:), 'b-', 'LineWidth', 1.6);
grid on; xlabel('Time [s]'); ylabel('Speed v [m/s]');
%title('Speed Tracking (Discrete-Time PID with Conditional Anti-Windup)');
legend('v_{ref}[k]','v[k]','Location','Best');

% Acceleration command
figure('Name','Acceleration Command','Color','w');
plot(t, a_cmd, 'r-', 'LineWidth', 1.6); grid on;
xlabel('Time [s]'); ylabel('a[k] [m/s^2]');
%title('Longitudinal Acceleration Command (Saturated)');

% Straight path (for completeness)
figure('Name','Path (Straight Line)','Color','w');
plot(x_log(1,:), x_log(2,:), 'b-', 'LineWidth', 1.6);
axis equal; grid on; xlabel('x [m]'); ylabel('y [m]');
%title('Vehicle Path (Straight Line, \delta_k = 0)');

% Speed error
figure('Name','Speed Error','Color','w');
plot(t, ev, 'LineWidth', 1.6); grid on;
xlabel('Time [s]'); ylabel('e_v[k] [m/s]');
%title('Speed Error Over Time');

% (Optional) Save figures to files as in the chapter
% saveas(1, 'figures/CH05/FigEx_PID_speed_straight_profile.jpg');
% saveas(2, 'figures/CH05/FigEx_PID_speed_straight_accel.jpg');
% saveas(3, 'figures/CH05/FigEx_PID_speed_straight_path.jpg');
% saveas(4, 'figures/CH05/FigEx_PID_speed_straight_error.jpg');
