% Example 7.1 (weight study): Path-indexed TV-LQR with a kinematic bicycle
% model, repeated for three weight sets to show the tracking/effort trade-off:
%   (1) baseline Q, R   (2) 10*Q (tighter tracking)   (3) 10*R (gentler inputs)
% Commands are limited to a in [-3, 2] m/s^2 and |delta| <= 30 deg (as in Chapter 6).
% Produces two figures for the book and prints RMS/peak metrics.
clc; clear; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

%% Exercise 7.3 option (false reproduces the constrained worked example)
remove_actuator_limits = false; % set true for the new end-of-chapter exercise

%% Common set-up (identical to Ex7_1_LQR_kinematic.m)
L = 2.8; Ts = 0.05; T_final = 200;
N = floor(T_final/Ts); Tspan = (0:N-1)*Ts;
[WP, velocities, headings] = waypoints_with_velocity_and_heading(N);
reference = [WP; unwrap(headings); velocities];
lim.a_min = -3; lim.a_max = 2;      % acceleration bounds [m/s^2]
lim.delta_max = deg2rad(30);        % steering bound [rad]

% Nominal inputs from path derivatives (waypoints are spaced by distance):
% a_nom = v dv/ds, delta_nom = atan(L*kappa), both kept within the limits
s_ref = [0, cumsum(vecnorm(diff(WP, 1, 2)))];
kappa = gradient(reference(3,:), s_ref);
dv_ds = gradient(reference(4,:), s_ref);
U_nom = zeros(2,N);
U_nom(1,:) = min(max(reference(4,:).*dv_ds, lim.a_min), lim.a_max);
U_nom(2,:) = min(max(atan(L*kappa), -lim.delta_max), lim.delta_max);

% Retain the same nominal trajectory; change applied-command limits only.
if remove_actuator_limits
    lim.a_min = -Inf; lim.a_max = Inf;
    lim.delta_max = Inf;
end

Q0 = diag([1, 1, 100, 10]);  R0 = diag([1, 1]);
cases = { 'Baseline Q, R',  Q0,    R0;
          '10Q',            10*Q0, R0;
          '10R',            Q0,    10*R0 };
styles = {'b-','r--','k-.'};

for c = 1:size(cases,1)
    res(c) = run_tvlqr(cases{c,2}, cases{c,3}, WP, headings, reference, U_nom, L, Ts, N, lim); %#ok<SAGROW>
end

%% Metrics
fprintf('\n%-14s %9s %9s %9s %9s %11s %9s %9s\n','Case','RMS CTE','max|CTE|','RMS e_th','RMS e_v','max|delta|','t_sat,a','lap time');
fprintf('%-14s %9s %9s %9s %9s %11s %9s %9s\n','','(m)','(m)','(deg)','(m/s)','(deg)','(s)','(s)');
for c = 1:numel(res)
    r = res(c);
    tsat = Ts*sum(r.u(1,:) >= lim.a_max - 1e-9 | r.u(1,:) <= lim.a_min + 1e-9);
    fprintf('%-14s %9.3f %9.3f %9.3f %9.3f %11.2f %9.2f %9.2f\n', cases{c,1}, ...
        rms(r.cte), max(abs(r.cte)), rms(rad2deg(r.the)), rms(r.ev), ...
        max(abs(rad2deg(r.u(2,:)))), tsat, r.t(end));
end

%% Plot defaults
set(groot,'defaultFigureColor','w','defaultAxesFontName','Helvetica', ...
    'defaultAxesFontSize',16,'defaultAxesFontWeight','bold', ...
    'defaultAxesLineWidth',1.5,'defaultLineLineWidth',2.5, ...
    'defaultLegendFontSize',13,'defaultLegendFontWeight','bold');
figdir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'figures', 'CH07');
if remove_actuator_limits, figdir = fullfile(figdir,'exercise7_3_unconstrained'); end
if ~exist(figdir, 'dir'), mkdir(figdir); end

%% Figure A: tracking errors
fA = figure('Position',[100 100 900 800]);
subplot(3,1,1); hold on;
for c = 1:numel(res), plot(res(c).t, res(c).cte, styles{c}); end
ylabel('CTE (m)'); grid on; box on; axis tight;
legend(cases(:,1),'Location','best');
subplot(3,1,2); hold on;
for c = 1:numel(res), plot(res(c).t, rad2deg(res(c).the), styles{c}); end
ylabel('Heading Error (deg)'); grid on; box on; axis tight;
subplot(3,1,3); hold on;
for c = 1:numel(res), plot(res(c).t, res(c).ev, styles{c}); end
ylabel('Speed Error (m/s)'); xlabel('Time (s)'); grid on; box on; axis tight;
exportgraphics(fA, fullfile(figdir,'FigEx7_1_errors.jpeg'), 'Resolution',300);

%% Figure B: control inputs
fB = figure('Position',[100 100 900 600]);
subplot(2,1,1); hold on;
for c = 1:numel(res), plot(res(c).t, res(c).u(1,:), styles{c}); end
ylabel('Acceleration (m/s^2)'); grid on; box on; axis tight;
legend(cases(:,1),'Location','best');
subplot(2,1,2); hold on;
for c = 1:numel(res), plot(res(c).t, rad2deg(res(c).u(2,:)), styles{c}); end
ylabel('Steering Angle (deg)'); xlabel('Time (s)'); grid on; box on; axis tight;
exportgraphics(fB, fullfile(figdir,'FigEx7_1_commands.jpeg'), 'Resolution',300);

%% ------------------------------------------------------------------------
function r = run_tvlqr(Q, R, WP, headings, reference, U_nom, L, Ts, N, lim)
X = zeros(4,N); X(:,1) = [WP(:,1); 0; 0];
U = zeros(2,N); cte_v = zeros(1,N); the_v = zeros(1,N); ev_v = zeros(1,N);
for k = 1:N-1
    [nxt, pas, i_k, ~] = nextWP(WP, X(1:2,k));
    if k == 1, first = nxt; end
    cte_v(k) = cte(nxt, pas, X(1:2,k));
    the_v(k) = mod(headings(i_k) - X(3,k) + pi, 2*pi) - pi;
    psi = reference(3,i_k); v_x = reference(4,i_k); delta = U_nom(2,i_k);
    A = [1 0 -Ts*v_x*sin(psi) Ts*cos(psi);
         0 1  Ts*v_x*cos(psi) Ts*sin(psi);
         0 0  1               Ts*tan(delta)/L;
         0 0  0               1];
    B = [0 0; 0 0; 0 (Ts/L)*v_x*sec(delta)^2; Ts 0];
    K = dlqr(A, B, Q, R);
    xt = [X(1:2,k) - reference(1:2,i_k);
          wrapToPi(X(3,k) - reference(3,i_k));
          X(4,k) - reference(4,i_k)];
    ev_v(k) = xt(4);
    U(:,k) = U_nom(:,i_k) - K*xt;
    U(1,k) = min(max(U(1,k), lim.a_min), lim.a_max);            % acceleration limit
    U(2,k) = min(max(U(2,k), -lim.delta_max), lim.delta_max);   % steering limit
    X(:,k+1) = X(:,k) + Ts*[X(4,k)*cos(X(3,k)); X(4,k)*sin(X(3,k));
                            X(4,k)/L*tan(U(2,k)); U(1,k)];
    if k > 100 && norm(first - X(1:2,k)) < 5, break; end
end
r.t = (0:k-1)*Ts; r.cte = cte_v(1:k); r.the = the_v(1:k);
r.ev = ev_v(1:k); r.u = U(:,1:k);
end
