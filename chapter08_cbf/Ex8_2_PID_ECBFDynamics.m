% Example 8.2: PID vs PID+ECBF Safety Filter (Dynamic Bicycle Model)
% Structure follows Ex6_4_PID_dynamics.m.
% A parallel PID+ECBF branch is added alongside the pure PID branch.
% Both share identical vehicle parameters, waypoints, and longitudinal control.
%
% The ECBF coefficients alpha and beta are obtained from the dynamic bicycle
% model of Chapter 4 (dynamic_bicycle_model.m), with two tires per axle:
%   e_ddot = v_t*vx*kappa
%            - (1/m)*[ (2*Cxr*Sr - CxA*vx^2)*sin(theta_e)
%                      + 2*Cxf*Sf*sin(theta_e + delta)
%                      + 2*Cyf*(delta - phi_f)*cos(theta_e + delta)
%                      - 2*Cyr*phi_r*cos(theta_e) ],
% where v_t = vx*cos(theta_e) - vy*sin(theta_e),
%       phi_f = (vy + Lf*omega)/vx,  phi_r = (vy - Lr*omega)/vx.

clc; clear; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

%% ─── VEHICLE PARAMETERS (identical to Ex6_4) ────────────────────────────
m   = 1500;
Lf  = 1.2;  Lr = 1.6;  L = Lf + Lr;
Iz  = m * ((0.5*(Lf+Lr))^2);
Cd  = 0.32;  A = 2.2;  rho = 1.225;
CxA = 0.5*Cd*A*rho;
Cxf = 60000;  Cxr = 120000;
Cyf = 30000;  Cyr = 30000;

%% ─── SIMULATION TIME (identical to Ex6_4) ───────────────────────────────
Ts      = 0.05;
T_final = 200;
N       = floor(T_final / Ts);
Tspan   = (0:N-1) * Ts;

% Integration of the vehicle model over one sampling period (input held):
%   true  -> fourth-order Runge-Kutta step
%   false -> forward Euler step
use_rk4 = true;

%% ─── REFERENCE TRAJECTORY (identical to Ex6_4) ─────────────────────────
[WP, velocities, headings] = waypoints_with_velocity_and_heading(N);
headings = unwrap(headings);
velocities = velocities*0.7;
reference = [WP; headings; velocities];
kappa_ref = reference_curvature(WP, headings);

%% ─── INITIAL STATES (identical to Ex6_4: zeros) ─────────────────────────
X_pid = zeros(6, N);   % PID branch:      [x y psi vx vy omega]
X_cbf = zeros(6, N);   % PID+ECBF branch: [x y psi vx vy omega]
X_pid(:,1) = [WP(:,1); headings(1); 0; 0; 0];
X_cbf(:,1) = [WP(:,1); headings(1); 0; 0; 0];

%% ─── LATERAL PID GAINS (identical to Ex6_4) ─────────────────────────────
Kp = 0.1;  Ki = 0;  Kd = 0.02;  Kt = 1;
int_pid = 0;   % lateral integral — PID branch
int_cbf = 0;   % lateral integral — ECBF branch

%% ─── LONGITUDINAL SPEED CONTROLLER ──────────────────────────────────────
% The waypoint speed profile contains abrupt changes.  A one-step command
% (v_ref - v_x)/Ts produces unrealistically large acceleration commands and
% amplifies the lateral-longitudinal coupling when the ECBF changes steering.
Kv_long = 0.1;      % first-order speed tracking gain (1/s)
ax_max  = 8.0;      % maximum commanded acceleration (m/s^2)
ax_min  = -10.0;    % maximum commanded braking acceleration (m/s^2)
Sr_max  = 0.20;     % rear slip-ratio saturation for numerical realism

%% ─── ECBF PARAMETERS ────────────────────────────────────────────────────
e_max   = 0.8;        % lateral safety bound (m)
lambda  = 5;
lambda0 = lambda^2;   % coefficient of h   in ECBF condition
lambda1 = 2*lambda;   % coefficient of hdot in ECBF condition
del_max = deg2rad(30);

%% ─── LOG ARRAYS ─────────────────────────────────────────────────────────
e_pid_log     = zeros(1, N);
e_cbf_log     = zeros(1, N);
theta_e_pid   = zeros(1, N);
theta_e_cbf   = zeros(1, N);
U_pid         = zeros(3, N);
U_cbf         = zeros(3, N);
h_log         = zeros(1, N);
d_delta_log   = zeros(1, N);
alpha_log     = zeros(1, N);
beta_log      = zeros(1, N);
v_ref_pid_log = zeros(1, N);
v_ref_cbf_log = zeros(1, N);
ecbf_active_log = false(1, N);
n_infeasible  = 0;     % number of steps at which the ECBF-QP has no solution

k_end = N - 1;

%% ─── MAIN LOOP ───────────────────────────────────────────────────────────
for k = 1:N-1

    % ══════════════════════════════════════════════════════════════════════
    %  PID BRANCH  (identical flow to Ex6_4)
    % ══════════════════════════════════════════════════════════════════════
    x_p   = X_pid(1,k);  y_p  = X_pid(2,k);
    psi_p = X_pid(3,k);  vx_p = X_pid(4,k);

    [next_p(:,k), pass_p(:,k), ni_p(k), ~] = nextWP(WP, [x_p; y_p]);

    e_pid_log(k)   = cte(next_p(:,k), pass_p(:,k), [x_p; y_p]);
    % Controller heading error uses the stabilizing PID convention
    % psi_ref - psi. The ECBF dynamics below use theta_e = psi - psi_ref,
    % matching the convention in the chapter text.
    theta_e_pid(k) = wrapToPiLocal(headings(ni_p(k)) - psi_p);

    int_pid = int_pid + e_pid_log(k)*Ts;
    if k == 1
        delta_p = Kp*e_pid_log(k) + Ki*int_pid + Kt*theta_e_pid(k);
    else
        delta_p = Kp*e_pid_log(k) + Ki*int_pid ...
                + Kd*(e_pid_log(k)-e_pid_log(k-1))/Ts + Kt*theta_e_pid(k);
    end
    delta_p = min(max(delta_p, -del_max), del_max);

    % Longitudinal speed tracking.  The acceleration command is bounded so
    % abrupt reference-speed changes do not dominate the lateral comparison.
    v_ref_pid = velocities(ni_p(k));
    v_ref_pid_log(k) = v_ref_pid;
    accel_p = Kv_long * (v_ref_pid - vx_p);
    accel_p = min(max(accel_p, ax_min), ax_max);
    Sr_p = (m * accel_p + CxA * vx_p^2) / (2*Cxr);
    Sr_p = min(max(Sr_p, -Sr_max), Sr_max);
    U_pid(:,k) = [delta_p; 0; Sr_p];

    X_pid(:,k+1) = plant_step(X_pid(:,k), U_pid(:,k), Ts, use_rk4, ...
               m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr);

    % ══════════════════════════════════════════════════════════════════════
    %  PID + ECBF BRANCH
    % ══════════════════════════════════════════════════════════════════════
    x_c   = X_cbf(1,k);  y_c   = X_cbf(2,k);
    psi_c = X_cbf(3,k);  vx_c  = X_cbf(4,k);
    vy_c  = X_cbf(5,k);  om_c  = X_cbf(6,k);
    vx_safe = max(abs(vx_c), 0.5);   % match dynamic_bicycle_model speed guard

    [next_c(:,k), pass_c(:,k), ni_c(k), ~] = nextWP(WP, [x_c; y_c]);

    e_cbf_log(k)   = cte(next_c(:,k), pass_c(:,k), [x_c; y_c]);
    % PID heading-error convention for the nominal controller.
    theta_e_cbf(k) = wrapToPiLocal(headings(ni_c(k)) - psi_c);

    % Nominal lateral PID (same structure as PID branch)
    int_cbf = int_cbf + e_cbf_log(k)*Ts;
    if k == 1
        delta_nom = Kp*e_cbf_log(k) + Ki*int_cbf + Kt*theta_e_cbf(k);
    else
        delta_nom = Kp*e_cbf_log(k) + Ki*int_cbf ...
                  + Kd*(e_cbf_log(k)-e_cbf_log(k-1))/Ts + Kt*theta_e_cbf(k);
    end
    delta_nom = min(max(delta_nom, -del_max), del_max);

    % Longitudinal command -- same bounded speed controller as PID branch.
    % It is computed before the ECBF terms because dvx/dt and dvy/dt enter
    % the full dynamic expression for e_ddot.
    v_ref_cbf = velocities(ni_c(k));
    v_ref_cbf_log(k) = v_ref_cbf;
    accel_c = Kv_long * (v_ref_cbf - vx_c);
    accel_c = min(max(accel_c, ax_min), ax_max);
    Sr_c = (m * accel_c + CxA * vx_c^2) / (2*Cxr);
    Sr_c = min(max(Sr_c, -Sr_max), Sr_max);
    Sf_c     = 0;

    % ── ECBF terms ────────────────────────────────────────────────────────
    % The chapter uses theta_e = psi - psi_ref and positive CTE on the
    % right-hand side of the path. The nominal PID above uses the opposite
    % heading-error sign only as a stabilizing controller convention.
    psi_ref   = headings(ni_c(k));
    kappa_c   = kappa_ref(ni_c(k));
    theta_dyn = wrapToPiLocal(psi_c - psi_ref);

    % Direct first derivative of the dynamic cross-track error:
    %   e_dot = -vx*sin(theta_e) - vy*cos(theta_e).
    e_dot = -vx_c*sin(theta_dyn) - vy_c*cos(theta_dyn);

    % alpha and beta of the affine approximation
    %   e_ddot ~= alpha*Delta_delta + beta,   Delta_delta = delta - delta_nom,
    % obtained from the dynamic bicycle model (two tires per axle).
    phi_f  = (vy_c + Lf*om_c)/vx_safe;       % front slip term
    phi_r  = (vy_c - Lr*om_c)/vx_safe;       % rear slip term
    s_e    = sin(theta_dyn);
    c_e    = cos(theta_dyn);
    v_t    = vx_c*c_e - vy_c*s_e;            % speed along the path tangent
    F_long = 2*Cxr*Sr_c - CxA*vx_c^2;        % rear traction minus drag (N)

    th_d   = theta_dyn + delta_nom;          % theta_e + delta_nom

    alpha = -(2/m) * ( (Cyf + Cxf*Sf_c)*cos(th_d) ...
                       - Cyf*(delta_nom - phi_f)*sin(th_d) );

    beta  = v_t*vx_c*kappa_c - (1/m) * ( F_long*s_e ...
              + 2*Cxf*Sf_c*sin(th_d) ...
              + 2*Cyf*(delta_nom - phi_f)*cos(th_d) ...
              - 2*Cyr*phi_r*c_e );

    alpha_log(k) = alpha;
    beta_log(k)  = beta;

    % Barrier value
    h_val    = e_max^2 - e_cbf_log(k)^2;
    h_log(k) = h_val;

    % ── ECBF QP ───────────────────────────────────────────────────────────
    % Solve for the steering correction Delta_delta:
    %   min  (Delta_delta)^2
    %   s.t. A_cbf * Delta_delta <= b_cbf
    %        -del_max <= delta_nom + Delta_delta <= del_max
    %
    % Constraint derivation:
    %   h_ddot + lambda1*h_dot + lambda0*h >= 0
    %   e_ddot ~= alpha*Delta_delta + beta
    %   -2(e_dot^2 + e*(alpha*Delta_delta+beta)) - 2*l1*e*e_dot + l0*h >= 0
    %   => 2*e*alpha*Delta_delta <= -2*(e_dot^2+e*beta) - 2*l1*e*e_dot + l0*h
    A_cbf = 2 * e_cbf_log(k) * alpha;

    b_cbf = -2*(e_dot^2 + e_cbf_log(k)*beta) ...
        - 2*lambda1*e_cbf_log(k)*e_dot ...
        + lambda0*h_val;

    lb = -del_max - delta_nom;
    ub =  del_max - delta_nom;

    [d_delta_safe, ~, flag] = quadprog(2, 0, A_cbf, b_cbf, ...
        [], [], lb, ub, [], optimoptions('quadprog','Display','off'));


    if flag == 1
        delta_c = delta_nom + d_delta_safe;
    else
        delta_c = delta_nom;              % fallback: nominal command
        n_infeasible = n_infeasible + 1;
    end

    d_delta_log(k) = delta_c - delta_nom;
    ecbf_active_log(k) = abs(d_delta_log(k)) > deg2rad(0.05);

    U_cbf(:,k) = [delta_c; Sf_c; Sr_c];

    X_cbf(:,k+1) = plant_step(X_cbf(:,k), U_cbf(:,k), Ts, use_rk4, ...
               m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr);

    % ── Termination (identical criterion to Ex6_4) ────────────────────────
    distanceToStart(k) = dis([next_c(1,1); next_c(2,1)], [x_c; y_c]);
    if distanceToStart(k) < 5 && k > 100
        fprintf('Stopping at step %d (t = %.1f s)\n', k, k*Ts);
        k_end = k;
        break;
    end
end

%% ─── PERFORMANCE METRICS ──────────────────────────────────────────────────
% Lateral tracking metrics for the completed portion of the simulation.
% A negative percentage change means the PID+ECBF branch reduced the metric
% relative to the nominal PID branch.
rmse_pid      = sqrt(mean(e_pid_log(1:k_end).^2));
rmse_ecbf     = sqrt(mean(e_cbf_log(1:k_end).^2));
max_abs_pid   = max(abs(e_pid_log(1:k_end)));
max_abs_ecbf  = max(abs(e_cbf_log(1:k_end)));

rmse_change_pct    = 100 * (rmse_ecbf - rmse_pid) / max(rmse_pid, eps);
max_abs_change_pct = 100 * (max_abs_ecbf - max_abs_pid) / max(max_abs_pid, eps);

fprintf('\n=== Performance Summary ===\n');
fprintf('Lateral RMSE  — PID: %.4f m   |   PID+ECBF: %.4f m   |   Change: %.1f%%\n', ...
        rmse_pid, rmse_ecbf, rmse_change_pct);
fprintf('Max |e|       — PID: %.4f m   |   PID+ECBF: %.4f m   |   Change: %.1f%%\n', ...
        max_abs_pid, max_abs_ecbf, max_abs_change_pct);
fprintf('Integration   — %s,  Ts = %.3f s,  lap time = %.2f s\n', ...
        ternary_str(use_rk4, 'RK4', 'forward Euler'), Ts, k_end*Ts);
fprintf('min h(e)      — PID+ECBF: %.4f\n', min(h_log(1:k_end)));
fprintf('ECBF active   — %.1f%% of the steps,  infeasible steps: %d\n', ...
        100*mean(ecbf_active_log(1:k_end)), n_infeasible);
fprintf('max |Delta_delta| = %.2f deg,  max |delta| (PID+ECBF) = %.2f deg,  max |delta| (PID) = %.2f deg\n', ...
        rad2deg(max(abs(d_delta_log(1:k_end)))), ...
        rad2deg(max(abs(U_cbf(1,1:k_end)))), rad2deg(max(abs(U_pid(1,1:k_end)))));
fprintf('Speed (PID+ECBF) — max: %.2f m/s,  desired on straights/arcs: %.1f / %.1f m/s\n', ...
        max(X_cbf(4,1:k_end)), max(velocities), min(velocities));
act_idx = find(ecbf_active_log(1:k_end));
if ~isempty(act_idx)
    fprintf('First ECBF activation at t = %.2f s\n', Tspan(act_idx(1)));
end

%% ─── SAFETY CORRIDOR ────────────────────────────────────────────────────
% Left-normal unit vector at each waypoint
nx = -sin(headings);   ny = cos(headings);
ux = WP(1,:) + e_max*nx;   uy = WP(2,:) + e_max*ny;
lx = WP(1,:) - e_max*nx;   ly = WP(2,:) - e_max*ny;

ecbf_color = 'green';      % dark orange for textbook visibility
active_color = [1.00 0.86 0.35];    % light amber activation highlight


%% ─── PLOT STYLE (same as Example 8.1) ───────────────────────────────────
set(groot,'defaultFigureColor','w', ...
          'defaultAxesFontName','Helvetica', ...
          'defaultAxesFontSize',16, ...
          'defaultAxesFontWeight','bold', ...
          'defaultAxesLineWidth',1.5, ...
          'defaultLineLineWidth',3, ...
          'defaultLegendFontSize',14, ...
          'defaultLegendFontWeight','bold');

%% ─── FIGURE 1: TRAJECTORY WITH CORRIDORS ────────────────────────────────
fig_path = figure;
main_ax = gca;
hold(main_ax, 'on'); grid(main_ax, 'on'); axis(main_ax, 'equal'); box(main_ax, 'on');

fill([ux, fliplr(lx)], [uy, fliplr(ly)], [0.85 0.92 1], ...
     'EdgeColor','none','FaceAlpha',0.4, 'DisplayName','Safety area');
plot(WP(1,:), WP(2,:), 'b-',  'LineWidth',2,   'DisplayName','Reference');
plot(ux, uy,            'k--', 'LineWidth',1.2, 'DisplayName', ...
     sprintf('Safety bound'));
plot(lx, ly,            'k--', 'LineWidth',1.2, 'HandleVisibility','off');
plot(X_pid(1,1:k_end), X_pid(2,1:k_end), 'r--', 'LineWidth',2, ...
     'DisplayName','PID');
plot(X_cbf(1,1:k_end), X_cbf(2,1:k_end), 'g-.', 'Color', ecbf_color, ...
     'LineWidth',2.4, ...
     'DisplayName','PID+ECBF');

xlabel('X Position (m)');  ylabel('Y Position (m)');
title('Trajectory — PID vs PID+ECBF');
% Legend below the axes, so that it does not cover the track.
legend(main_ax, 'Location','southoutside', 'Orientation','horizontal');

% Inset: entry to the second arc (upper left), where the cross-track error
% of the PID controller is largest. The window is 13 m x 7 m, so the
% 1.6 m wide corridor is clearly resolved.
zoom_xlim = [-4, 9];
zoom_ylim = [101.3, 108.3];

rectangle(main_ax, 'Position', ...
    [zoom_xlim(1), zoom_ylim(1), diff(zoom_xlim), diff(zoom_ylim)], ...
    'EdgeColor', 'k', 'LineWidth', 1.0, 'LineStyle', '-', ...
    'HandleVisibility', 'off');

% The inset is placed inside the loop. Its rectangle is given in data units
% of the main axes, with the same aspect ratio as the zoom window, and is
% converted to figure units, so that it stays clear of the track for any
% figure size and legend placement.
inset_w    = 130;
inset_rect = [8, 20, inset_w, inset_w*diff(zoom_ylim)/diff(zoom_xlim)];

% Connect the zoom window to the inset.
plot(main_ax, [zoom_xlim(1), inset_rect(1)], ...
     [zoom_ylim(1), inset_rect(2) + inset_rect(4)], 'k:', ...
     'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(main_ax, [zoom_xlim(2), inset_rect(1) + inset_rect(3)], ...
     [zoom_ylim(1), inset_rect(2) + inset_rect(4)], 'k:', ...
     'LineWidth', 1.0, 'HandleVisibility', 'off');

[fx1, fy1] = dataToFigure(main_ax, inset_rect(1), inset_rect(2));
[fx2, fy2] = dataToFigure(main_ax, inset_rect(1) + inset_rect(3), ...
                                   inset_rect(2) + inset_rect(4));
inset_ax = axes('Parent', fig_path, 'Units', 'normalized', ...
    'Position', [fx1, fy1, fx2 - fx1, fy2 - fy1], 'Color', 'w');
hold(inset_ax, 'on'); grid(inset_ax, 'on'); axis(inset_ax, 'equal');
box(inset_ax, 'on');
fill(inset_ax, [ux, fliplr(lx)], [uy, fliplr(ly)], [0.85 0.92 1], ...
     'EdgeColor','none','FaceAlpha',0.4);
plot(inset_ax, WP(1,:), WP(2,:), 'b-', 'LineWidth',2);
plot(inset_ax, ux, uy, 'k--', 'LineWidth',1.2);
plot(inset_ax, lx, ly, 'k--', 'LineWidth',1.2);
plot(inset_ax, X_pid(1,1:k_end), X_pid(2,1:k_end), 'r--', ...
     'LineWidth',2);
plot(inset_ax, X_cbf(1,1:k_end), X_cbf(2,1:k_end), '-.', ...
     'Color', ecbf_color, 'LineWidth',2.4);
xlim(inset_ax, zoom_xlim); ylim(inset_ax, zoom_ylim);
%title(inset_ax, 'Zoomed view', 'FontSize', 8);
set(inset_ax, 'FontSize', 11, 'LineWidth', 1.2);
uistack(inset_ax, 'top');
axes(inset_ax);

%% ─── FIGURE 2: CTE AND BARRIER FUNCTION ────────────────────────────────
fig_cte = figure;
subplot(2,1,1); hold on; grid on; box on;
plot(Tspan(1:k_end), e_pid_log(1:k_end), 'b', 'LineWidth',1.5, ...
     'DisplayName','PID');
plot(Tspan(1:k_end), e_cbf_log(1:k_end), 'r--', ...
     'LineWidth',1.8, ...
     'DisplayName','PID+ECBF');
yline( e_max, 'k--', 'LineWidth',1.2, 'HandleVisibility','off');
yline(-e_max, 'k--', 'LineWidth',1.2, 'HandleVisibility','off');
ylabel('CTE (m)');
title('Cross-Track Error');
shadeActivation(Tspan(1:k_end), ecbf_active_log(1:k_end), active_color);
% Horizontal legend centered at (t, e) = (17.5 s, -0.4 m), which lies in the
% free band between the error curves and the lower bound.
lgd = legend('Orientation','horizontal');
[cx, cy] = dataToFigure(gca, 17.5, -0.4);
lgd.Units = 'normalized';
lp = lgd.Position;
lgd.Position = [cx - lp(3)/2, cy - lp(4)/2, lp(3), lp(4)];

subplot(2,1,2); hold on; grid on; box on;
plot(Tspan(1:k_end), h_log(1:k_end), 'b', ...
     'LineWidth',1.8);
yline(0, 'k--', 'LineWidth',1.2);
ylabel('h(e) = e_{max}^2 - e^2');
xlabel('Time (s)');
title('Barrier Function h(e)  [PID+ECBF — must stay \geq 0]');
shadeActivation(Tspan(1:k_end), ecbf_active_log(1:k_end), active_color);

% %% ─── FIGURE 3: STEERING ANGLE AND ECBF CORRECTION ───────────────────────
% %figure('Color','w','Position',[80 60 900 600]);
% figure;
% subplot(2,1,1); hold on; grid on; box on;
% plot(Tspan(1:k_end), rad2deg(U_pid(1,1:k_end)), 'r--', 'LineWidth',1.5, ...
%      'DisplayName','PID');
% plot(Tspan(1:k_end), rad2deg(U_cbf(1,1:k_end)), '-.', ...
%      'Color', ecbf_color, 'LineWidth',1.5, ...
%      'DisplayName','PID+ECBF');
% yline( rad2deg(del_max), 'k--', 'LineWidth',1, 'HandleVisibility','off');
% yline(-rad2deg(del_max), 'k--', 'LineWidth',1, 'HandleVisibility','off');
% ylabel('\delta (deg)');
% title('Steering Angle \delta');
% legend('Location','best');
% 
% subplot(2,1,2); hold on; grid on; box on;
% plot(Tspan(1:k_end), rad2deg(d_delta_log(1:k_end)), 'LineWidth',1.5);
% yline(0, 'k--');
% ylabel('\Delta\delta^* (deg)');
% xlabel('Time (s)');
% title('ECBF Steering Correction \Delta\delta^*  (non-zero = filter active)');

% %% ─── FIGURE 4: VELOCITIES ───────────────────────
fig_speed = figure; hold on;
plot(Tspan(1:k_end), v_ref_pid_log(1:k_end), 'k:', 'LineWidth', 2.0, ...
    'DisplayName', 'Desired');
plot(Tspan(1:k_end), X_pid(4,1:k_end), 'b-', 'LineWidth', 1.5, ...
    'DisplayName', 'PID');
plot(Tspan(1:k_end), X_cbf(4,1:k_end), 'r--', 'LineWidth', 1.5, ...
    'DisplayName', 'PID+ECBF');
xlim('tight');
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)');
legend('Location','southeast'); grid on; box on;

%% ─── SAVE THE THREE FIGURES USED IN THE BOOK (figures/CH08) ─────────────
% Saved only for the Runge-Kutta run, which is the one reported in the book.
if use_rk4
    figdir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'figures', 'CH08');
    exportgraphics(fig_path,  fullfile(figdir, 'FigEx7_2_path.jpeg'),  'Resolution', 300);
    exportgraphics(fig_cte,   fullfile(figdir, 'FigEx7_2_cte.jpeg'),   'Resolution', 300);
    exportgraphics(fig_speed, fullfile(figdir, 'FigEx7_2_speed.jpeg'), 'Resolution', 300);
end

function X_next = plant_step(X, U, Ts, use_rk4, m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr)
%PLANT_STEP Advance the dynamic bicycle model by one sampling period with
%the input U held constant.
    f = @(x) dynamic_bicycle_model(0, x, U, m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr);
    if use_rk4
        k1 = f(X);
        k2 = f(X + 0.5*Ts*k1);
        k3 = f(X + 0.5*Ts*k2);
        k4 = f(X + Ts*k3);
        X_next = X + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4);
    else
        X_next = X + Ts*f(X);
    end
end

function [fx, fy] = dataToFigure(ax, x, y)
%DATATOFIGURE Normalized figure coordinates of the data point (x, y) of ax.
%With "axis equal" the plot box is centered inside the axes position.
    drawnow;
    fig = ancestor(ax, 'figure');
    u = ax.Units;   ax.Units  = 'pixels';  ap = ax.Position;   ax.Units  = u;
    u = fig.Units;  fig.Units = 'pixels';  fp = fig.Position;  fig.Units = u;
    xl = xlim(ax);  yl = ylim(ax);
    sx = ap(3)/diff(xl);                 % pixels per data unit
    sy = ap(4)/diff(yl);
    if strcmp(ax.DataAspectRatioMode, 'manual')      % axis equal
        s  = min(sx, sy);
        ap(1) = ap(1) + (ap(3) - s*diff(xl))/2;
        ap(2) = ap(2) + (ap(4) - s*diff(yl))/2;
        sx = s;  sy = s;
    end
    fx = (ap(1) + sx*(x - xl(1)))/fp(3);
    fy = (ap(2) + sy*(y - yl(1)))/fp(4);
end

function s = ternary_str(cond, a, b)
    if cond, s = a; else, s = b; end
end

function angle = wrapToPiLocal(angle)
%WRAPTOPILOCAL Normalize angle to [-pi, pi] without toolbox dependencies.
    angle = mod(angle + pi, 2*pi) - pi;
end

function shadeActivation(t, active, color)
%SHADEACTIVATION Add pale background bands where the ECBF filter is active.
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

function kappa = reference_curvature(WP, headings)
%REFERENCE_CURVATURE Signed path curvature d(psi_ref)/ds.
    ds = sqrt(sum(diff(WP, 1, 2).^2, 1));
    s = [0, cumsum(max(ds, eps))];
    kappa = gradient(headings, s);
    kappa(~isfinite(kappa)) = 0;
end
