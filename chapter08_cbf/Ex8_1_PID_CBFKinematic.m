% Example 8.1: PID + Kinematic CBF Safety Filter for Lateral Control
%              with Longitudinal PID Speed Control
%
% Runs two parallel simulations:
%   (1) Pure PID     - lateral PID steering + longitudinal PID speed (baseline)
%   (2) PID + CBF    - same PID nominal command filtered by the kinematic CBF-QP:
%
%       min   (delta - delta_nom)^2
%       s.t.  2*e*vx*delta + lambda*(e_max^2 - e^2) >= 0   [CBF condition]
%             -del_max <= delta <= del_max                  [actuator limit]
%
% Both vehicles start with a 0.35 m lateral offset so the CBF filter
% is demonstrably active at the start of the simulation.

clear; clc; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

%% Parameters
L       = 1.08;            % Wheelbase (m)
N       = 1000;            % Simulation steps
del_max = deg2rad(30);     % Max steering angle (rad)
Ts      = 0.05;            % Sampling time (s)

%% Reference trajectory
n = 100;
[WP, velocities, headings] = waypoints_with_velocity_and_heading(n);
reference = [WP; unwrap(headings); velocities];

%% Lateral PID gains (CTE + heading)
Kp = 0.1;    % CTE proportional
Ki = 0.0;    % CTE integral
Kd = 0.0;    % CTE derivative
Kt = 1.0;    % Heading proportional

%% Longitudinal PID gains (speed control)
KvP   =  0.8;    % Speed proportional
KvI   =  0.2;    % Speed integral
KvD   =  0.05;   % Speed derivative
a_min = -3.0;    % Max deceleration (m/s^2)
a_max =  2.0;    % Max acceleration (m/s^2)

%% CBF parameters (lateral safety)
e_max  = 0.4;    % Maximum allowable cross-track error (m)
lambda = 1.0;    % CBF decay rate lambda > 0 (higher => more conservative)

%% Initial states — 0.35 m lateral offset to activate CBF at start
x_pid = [WP(1,1); WP(2,1) + 0.35; 0; 0];   % [X; Y; psi; vx]
x_cbf = [WP(1,1); WP(2,1) + 0.35; 0; 0];

%% Control input arrays  [row 1: acceleration;  row 2: steering angle]
u_pid = [zeros(1, N);  zeros(1, N)];
u_cbf = [zeros(1, N);  zeros(1, N)];

%% Lateral integrators
int_cte_pid = 0;
int_cte_cbf = 0;

%% Longitudinal integrators and previous speed errors
iv_pid = 0;   ev_prev_pid = 0;
iv_cbf = 0;   ev_prev_cbf = 0;

%% Simulation loop
for k = 1:N

    % ── Nearest waypoints ────────────────────────────────────────────────
    [next_pid(:,k), pass_pid(:,k), ni_pid(k), ~] = nextWP(WP, [x_pid(1,k); x_pid(2,k)]);
    [next_cbf(:,k), pass_cbf(:,k), ni_cbf(k), ~] = nextWP(WP, [x_cbf(1,k); x_cbf(2,k)]);

    % ── Desired states ───────────────────────────────────────────────────
    x_desired_pid(3,k) = reference(3, ni_pid(k));   % desired heading
    x_desired_pid(4,k) = reference(4, ni_pid(k));   % desired speed

    x_desired_cbf(3,k) = reference(3, ni_cbf(k));
    x_desired_cbf(4,k) = reference(4, ni_cbf(k));

    % ── Lateral errors ───────────────────────────────────────────────────
    e_cte_pid(k)  = cte(next_pid(:,k), pass_pid(:,k), [x_pid(1,k); x_pid(2,k)]);
    theta_e_pid(k)= atan2(sin(x_desired_pid(3,k) - x_pid(3,k)), ...
                          cos(x_desired_pid(3,k) - x_pid(3,k)));

    e_cte_cbf(k)  = cte(next_cbf(:,k), pass_cbf(:,k), [x_cbf(1,k); x_cbf(2,k)]);
    theta_e_cbf(k)= atan2(sin(x_desired_cbf(3,k) - x_cbf(3,k)), ...
                          cos(x_desired_cbf(3,k) - x_cbf(3,k)));

    % ── Longitudinal PID — PID branch ────────────────────────────────────
    ev_pid(k) = x_desired_pid(4,k) - x_pid(4,k);
    dev_pid   = (k==1)*0 + (k>1)*(ev_pid(k) - ev_prev_pid)/Ts;
    iv_try    = iv_pid + ev_pid(k)*Ts;
    a_unsat   = KvP*ev_pid(k) + KvI*iv_try + KvD*dev_pid;
    a_sat     = min(max(a_unsat, a_min), a_max);
    if abs(a_unsat - a_sat) < 1e-9,  iv_pid = iv_try;  end   % anti-windup
    u_pid(1,k)  = a_sat;
    ev_prev_pid = ev_pid(k);

    % ── Longitudinal PID — CBF branch ────────────────────────────────────
    ev_cbf(k) = x_desired_cbf(4,k) - x_cbf(4,k);
    dev_cbf   = (k==1)*0 + (k>1)*(ev_cbf(k) - ev_prev_cbf)/Ts;
    iv_try    = iv_cbf + ev_cbf(k)*Ts;
    a_unsat   = KvP*ev_cbf(k) + KvI*iv_try + KvD*dev_cbf;
    a_sat     = min(max(a_unsat, a_min), a_max);
    if abs(a_unsat - a_sat) < 1e-9,  iv_cbf = iv_try;  end   % anti-windup
    u_cbf(1,k)  = a_sat;
    ev_prev_cbf = ev_cbf(k);

    % ── Lateral PID — PID branch (pure PID steering) ─────────────────────
    if k == 1
        u_pid(2,k) = Kp*e_cte_pid(k) + Ki*int_cte_pid + Kt*theta_e_pid(k);
    else
        int_cte_pid = int_cte_pid + e_cte_pid(k)*Ts;
        u_pid(2,k)  = Kp*e_cte_pid(k) + Ki*int_cte_pid ...
                    + Kd*(e_cte_pid(k) - e_cte_pid(k-1))/Ts + Kt*theta_e_pid(k);
    end
    u_pid(2,k) = max(min(u_pid(2,k), del_max), -del_max);

    % ── Lateral PID — CBF branch nominal command ──────────────────────────
    if k == 1
        delta_nom(k) = Kp*e_cte_cbf(k) + Ki*int_cte_cbf + Kt*theta_e_cbf(k);
    else
        int_cte_cbf  = int_cte_cbf + e_cte_cbf(k)*Ts;
        delta_nom(k) = Kp*e_cte_cbf(k) + Ki*int_cte_cbf ...
                     + Kd*(e_cte_cbf(k) - e_cte_cbf(k-1))/Ts + Kt*theta_e_cbf(k);
    end

    % ── Kinematic CBF-QP ─────────────────────────────────────────────────
    %  Barrier:  h(e) = e_max^2 - e^2       (safe set: |e| <= e_max)
    %  CBF cond: hdot + lambda*h >= 0
    %            2*e*vx*delta + lambda*(e_max^2 - e^2) >= 0
    %  quadprog  A_c*delta <= b_c   with
    %            A_c = -2*e*vx,   b_c = lambda*(e_max^2 - e^2)
    h(k)  = e_max^2 - e_cte_cbf(k)^2;
    A_c   = -2 * e_cte_cbf(k) * x_cbf(4,k);
    b_c   =  lambda * h(k);

    [delta_safe, ~, flag] = quadprog(2, -2*delta_nom(k), A_c, b_c, ...
                                     [], [], -del_max, del_max, [], ...
                                     optimoptions('quadprog','Display','off'));
    if flag == 1
        u_cbf(2,k) = delta_safe;
    else
        u_cbf(2,k) = max(min(delta_nom(k), del_max), -del_max);  % fallback
    end


    % ── State update (kinematic bicycle model) ───────────────────────────
    x_pid(:, k+1) = x_pid(:, k) + Ts * [x_pid(4,k) * cos(x_pid(3,k));
                                          x_pid(4,k) * sin(x_pid(3,k));
                                          (x_pid(4,k) / L) * tan(u_pid(2,k));
                                          u_pid(1,k)];

    x_cbf(:, k+1) = x_cbf(:, k) + Ts * [x_cbf(4,k) * cos(x_cbf(3,k));
                                          x_cbf(4,k) * sin(x_cbf(3,k));
                                          (x_cbf(4,k) / L) * tan(u_cbf(2,k));
                                          u_cbf(1,k)];

    % ── Termination ──────────────────────────────────────────────────────
    distanceToStart(k) = dis([WP(1,1); WP(2,1)], [x_cbf(1,k); x_cbf(2,k)]);
    if distanceToStart(k) < 5 && k > 100
        fprintf('Stopping simulation at step %d because distance is less than the threshold.\n', k);
        break;
    end

    % Live plot (optional — comment out to speed up simulation)
    % plot(x_pid(1,k), x_pid(2,k), 'r*'); hold on;
    % plot(x_cbf(1,k), x_cbf(2,k), 'g*'); drawnow;
end

%% Performance metrics
t = (0:Ts:(k-1)*Ts);
fprintf('\n=== Performance Summary ===\n');
fprintf('Lateral RMSE  — PID: %.4f m   |   PID+CBF: %.4f m\n', ...
        sqrt(mean(e_cte_pid(1:k).^2)), sqrt(mean(e_cte_cbf(1:k).^2)));
fprintf('Max |CTE|     — PID: %.4f m   |   PID+CBF: %.4f m\n', ...
        max(abs(e_cte_pid(1:k))), max(abs(e_cte_cbf(1:k))));
fprintf('Speed  RMSE   — PID: %.4f m/s |   PID+CBF: %.4f m/s\n', ...
        sqrt(mean(ev_pid(1:k).^2)), sqrt(mean(ev_cbf(1:k).^2)));

%% Set Global defaults
set(groot,'defaultFigureColor','w', ...
          'defaultAxesFontName','Helvetica', ...
          'defaultAxesFontSize',16, ...
          'defaultAxesFontWeight','bold', ...
          'defaultAxesLineWidth',1.5, ...
          'defaultLineLineWidth',3, ...
          'defaultLegendFontSize',14, ...
          'defaultLegendFontWeight','bold');

%% Figure 1 — Trajectory comparison with magnified inset
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [100, 100, 900, 600]);
nx = -sin(headings);   ny = cos(headings);

% ── Main axes ────────────────────────────────────────────────────────────────
main_ax = axes('Position', [0.10, 0.11, 0.85, 0.82]);
hold(main_ax, 'on'); grid(main_ax, 'on');

plot(main_ax, reference(1,:), reference(2,:), 'b-',  'LineWidth', 2, 'DisplayName', 'Reference');
plot(main_ax, WP(1,:) + e_max*nx, WP(2,:) + e_max*ny, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Safety Bound');
plot(main_ax, WP(1,:) - e_max*nx, WP(2,:) - e_max*ny, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot(main_ax, x_pid(1,1:k), x_pid(2,1:k), 'r--', 'LineWidth', 2, 'DisplayName', 'PID');
plot(main_ax, x_cbf(1,1:k), x_cbf(2,1:k), 'g-.', 'LineWidth', 2, 'DisplayName', 'PID + CBF');

xlabel(main_ax, 'X Position (m)');
ylabel(main_ax, 'Y Position (m)');
legend(main_ax, 'Location', 'best');
axis(main_ax, 'equal'); box(main_ax, 'on');

% ── Define zoom region (YOU CAN TUNE THIS) ───────────────────────────────────
x_zoom = [151 155];   % <-- adjust
y_zoom = [14 18];    % <-- adjust

% Draw rectangle on main plot
rectangle(main_ax, ...
    'Position', [x_zoom(1), y_zoom(1), diff(x_zoom), diff(y_zoom)], ...
    'EdgeColor', 'k', 'LineStyle', ':', 'LineWidth', 1.5);

% ── Inset axes ───────────────────────────────────────────────────────────────
inset_ax = axes('Position', [0.35, 0.35, 0.30, 0.30]);
hold(inset_ax, 'on'); grid(inset_ax, 'on');
box(inset_ax, 'on');

% Replot same data
plot(inset_ax, reference(1,:), reference(2,:), 'b-',  'LineWidth', 2);
plot(inset_ax, WP(1,:) + e_max*nx, WP(2,:) + e_max*ny, 'k--', 'LineWidth', 1.2);
plot(inset_ax, WP(1,:) - e_max*nx, WP(2,:) - e_max*ny, 'k--', 'LineWidth', 1.2);
plot(inset_ax, x_pid(1,1:k), x_pid(2,1:k), 'r--', 'LineWidth', 2);
plot(inset_ax, x_cbf(1,1:k), x_cbf(2,1:k), 'g-.', 'LineWidth', 2);

% Zoom limits
xlim(inset_ax, x_zoom);
ylim(inset_ax, y_zoom);
%axis(inset_ax, 'equal');

% % ── Add annotations inside inset ─────────────────────────────────────────────
% text(inset_ax, mean(x_zoom), y_zoom(2), ...
%     'Zoomed Region', 'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
% 
% % Optional: highlight difference between PID and CBF
% text(inset_ax, x_zoom(1)+0.5, y_zoom(1)+0.5, ...
%     'CBF stays within bounds', ...
%     'Color', 'g', 'FontSize', 9);
% 
% text(inset_ax, x_zoom(1)+0.5, y_zoom(1)+1.5, ...
%     'PID may violate constraint', ...
%     'Color', 'r', 'FontSize', 9);

% ── Optional: connect inset to main plot ─────────────────────────────────────
annotation('line', [0.65 0.78], [0.52 0.40], 'LineStyle', '--');
annotation('line', [0.55 0.45], [0.55 0.40], 'LineStyle', '--');

% %% Figure 2 — States: heading and velocity
% figure(2);
% subplot(2,1,1); hold on;
% plot(t, rad2deg(x_desired_pid(3,1:k)), 'b-',  'LineWidth', 2);
% plot(t, rad2deg(x_pid(3,1:k)),         'r--', 'LineWidth', 2);
% plot(t, rad2deg(x_cbf(3,1:k)),         'g-.', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel('Heading (deg)');
% legend('Reference','PID','PID + CBF','Location','best');
% grid on; box on; axis tight;
% 
% subplot(2,1,2); hold on;
% plot(t, x_desired_pid(4,1:k), 'b-',  'LineWidth', 2);
% plot(t, x_pid(4,1:k),         'r--', 'LineWidth', 2);
% plot(t, x_cbf(4,1:k),         'g-.', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel('Velocity (m/s)');
% legend('Reference','PID','PID + CBF','Location','best');
% grid on; box on; axis tight;

%% Figure 3 — Errors: CTE, heading, speed, and barrier function
figure(3); %subplot(2,1,1); 
hold on;
plot(t, e_cte_pid(1:k), 'r--', 'LineWidth', 2);
plot(t, e_cte_cbf(1:k), 'g-.', 'LineWidth', 2);
yline( e_max, 'k--', 'LineWidth', 1.5);
yline(-e_max, 'k--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('CTE (m)');
legend('PID','PID + CBF','Safety Bound','Location','best');
grid on; box on; axis tight;

% subplot(2,1,2); hold on;
% plot(t, rad2deg(theta_e_pid(1:k)), 'r--', 'LineWidth', 2);
% plot(t, rad2deg(theta_e_cbf(1:k)), 'g-.', 'LineWidth', 2);
% ylabel('Heading Error (deg)');
% legend('PID','PID + CBF','Location','best');
% grid on; box on; axis tight;

% subplot(3,1,3); hold on;
% plot(t, ev_pid(1:k), 'r--', 'LineWidth', 2);
% plot(t, ev_cbf(1:k), 'g-.', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel('Speed Error (m/s)');
% legend('PID','PID + CBF','Location','best');
% grid on; box on; axis tight;

% %% Figure 4 — Control inputs: steering and acceleration
% figure(4);
% subplot(2,1,1); hold on;
% plot(t, rad2deg(delta_nom(1:k)),    'b--', 'LineWidth', 2);
% plot(t, rad2deg(u_pid(2,1:k)),      'r-',  'LineWidth', 2);
% plot(t, rad2deg(u_cbf(2,1:k)),      'g-.', 'LineWidth', 2);
% yline( rad2deg(del_max), 'k--', 'LineWidth', 1);
% yline(-rad2deg(del_max), 'k--', 'LineWidth', 1);
% ylabel('Steering Angle (deg)');
% legend('Nominal \delta_{nom}','PID','Safe \delta^* (PID+CBF)','Location','best');
% grid on; box on; axis tight;
% 
% subplot(2,1,2); hold on;
% plot(t, u_pid(1,1:k), 'r-',  'LineWidth', 2);
% plot(t, u_cbf(1,1:k), 'g-.', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
% legend('PID','PID + CBF','Location','best');
% grid on; box on; axis tight;

% %% ─── FIGURE 4: VELOCITIES ───────────────────────
figure; hold on;
plot(t, x_desired_pid(4,1:k), 'b-',  'LineWidth', 2);
plot(t, x_pid(4,1:k),         'r--', 'LineWidth', 2);
plot(t, x_cbf(4,1:k),         'g-.', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Reference','PID','PID + CBF','Location','best');
xlim('tight');
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)');
legend('Location','best'); grid on; box on;

%% Figure 5 — CBF analysis: barrier function value and steering correction
active_color = [1.00 0.86 0.35];   % light amber activation highlight
correction = rad2deg(u_cbf(2,1:k) - min(max(delta_nom(1:k), -del_max), del_max));
cbf_active = abs(correction) > 0.05; % nonzero CBF steering intervention (deg)

figure(5);
subplot(2,1,1); hold on;
plot(t, h(1:k), 'b-', 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('h(e)');
legend('Barrier function  h(e) = e_{max}^2 - e^2','h = 0  (boundary)','Location','best');
grid on; box on; axis tight;
shadeActivation(t, cbf_active, active_color);

subplot(2,1,2); hold on;
plot(t, correction, 'k-', 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('\Delta\delta^{*}  (deg)');
legend('CBF correction \Delta\delta^*','Location','best');
grid on; box on; axis tight;
shadeActivation(t, cbf_active, active_color);

function shadeActivation(t, active, color)
%SHADEACTIVATION Add pale background bands where the safety filter is active.
    active = active(:)' ~= 0;
    if ~any(active)
        return;
    end

    yl = ylim;
    edges = diff([false, active, false]);
    starts = find(edges == 1);
    stops = find(edges == -1) - 1;

    for i = 1:numel(starts)
        x1 = t(starts(i));
        x2 = t(stops(i));
        if x2 <= x1
            x2 = x1 + eps;
        end
        p = patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], color, ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        uistack(p, 'bottom');
    end
    ylim(yl);
end
