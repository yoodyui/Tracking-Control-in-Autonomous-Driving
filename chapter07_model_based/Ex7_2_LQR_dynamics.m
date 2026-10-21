% Example 7.2: Path-Indexed TV-LQR with a Dynamic Bicycle Model
%Lateral + Longitudinal control
%Picking the closest WP in the loop
%Work very well for all velocity!!
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

%% Actuator limits (same acceleration range as in Chapter 6 and Example 7.1)
a_min = -3; a_max = 2;                 % longitudinal acceleration bounds [m/s^2]
delta_max = deg2rad(30);               % steering bound [rad]
Sr_min = m*a_min/(2*Cxr);              % rear slip ratio that gives a_min
Sr_max = m*a_max/(2*Cxr);              % rear slip ratio that gives a_max

%% Nominal inputs from the reference (path derivatives)
% The waypoints are spaced by distance, not by time: d(.)/dt = v * d(.)/ds
s_ref = [0, cumsum(vecnorm(diff(WP, 1, 2)))];        % arc length [m]
kappa = gradient(unwrap(headings), s_ref);           % curvature d(psi)/ds

%% Feasible nominal speed (forward-backward pass along the path)
% The curvature-based speed profile changes in steps, which no vehicle can
% follow. The nominal speed is therefore limited so that it can be reached
% with 80% of the acceleration range; the rest is left to the feedback.
margin = 0.8;
ds_c  = [diff(s_ref), norm(WP(:,1) - WP(:,N))];      % distance from waypoint i to i+1 (cyclic)
v_nom = velocities;                                  % start from the speed limit
for i = [1:N, 1:N]                                   % forward pass: acceleration limit
    j = mod(i, N) + 1;                               % (two laps make the profile periodic)
    v_nom(j) = min(v_nom(j), sqrt(v_nom(i)^2 + 2*margin*a_max*ds_c(i)));
end
for i = [N:-1:1, N:-1:1]                             % backward pass: braking limit
    j = mod(i, N) + 1;
    v_nom(i) = min(v_nom(i), sqrt(v_nom(j)^2 + 2*margin*abs(a_min)*ds_c(i)));
end
a_nom = v_nom .* gradient(v_nom, s_ref);             % a = v dv/ds
a_nom = min(max(a_nom, a_min), a_max);
Sf = 0;                                              % rear-wheel drive
U_nom = zeros(3,N);                                  % [delta; Sf; Sr]
U_nom(1,:) = min(max(atan(L*kappa), -delta_max), delta_max);
U_nom(2,:) = Sf;
U_nom(3,:) = m*a_nom/(2*Cxr);

%% Nominal states
[x_desired, y_desired] = deal(WP(1,:), WP(2,:));
psi_desired = unwrap(headings);
vx_desired = v_nom;                        % feasible nominal speed
vy_desired = zeros(1,N); % assuming nominal lateral velocity is zero
omega_desired = vx_desired .* kappa;       % nominal yaw rate = v * curvature
X_nom = [x_desired; y_desired; psi_desired; vx_desired; vy_desired; omega_desired];

%% Initialize states
X_lqr = zeros(6, N);
X_lqr(:,1) = [WP(:,1); psi_desired(1); vx_desired(1); 0; 0];   % start on the path at the nominal speed

%% LQR Weight Matrices
R = diag([1, 1, 0.1]);
Q = diag([10, 10, 500, 100, 1e-3, 100]);

%% Simulation Loop
for k = 1:N-1
   

    % Compute lateral error from LQR
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
    v_y = X_nom(5,next_index_lqr(k));
    omega = X_nom(6,next_index_lqr(k));

    %control inputs
    delta = U_nom(1,next_index_lqr(k));

   
    %% LQR design
    % Recompute Transformation Coefficients based on updated nominal parameters
    alpha13 = -(v_x * sin(psi) + v_y * cos(psi));
    alpha14 = cos(psi);
    alpha15 = -sin(psi);
    alpha23 = v_x * cos(psi) - v_y * sin(psi);
    alpha24 = sin(psi);
    alpha25 = cos(psi);
    alpha36 = 1;

    alpha44 = (1/m) * (-2 * Cyf * sin(delta) * (v_y + Lf * omega) / v_x^2 - 2 * CxA * v_x);
    alpha45 = (1/m) * (m * omega + 2 * Cyf * sin(delta) / v_x);
    alpha46 = (1/m) * (m * v_y + 2 * Cyf * sin(delta) * Lf / v_x);
    beta41 = (1/m) * (-2 * Cxf * Sf * sin(delta) - 2 * Cyf * (sin(delta) + (delta - (v_y + Lf * omega) / v_x) * cos(delta)));
    beta42 = (1/m) * (2 * Cxf * cos(delta));
    beta43 = (1/m) * (2 * Cxr);

    alpha54 = (1/m) * (-m * omega + 2 * Cyf * cos(delta) * (v_y + Lf * omega) / v_x^2 + 2 * Cyr * (v_y - Lr * omega) / v_x^2);
    alpha55 = (1/m) * (-2 * Cyf * cos(delta) / v_x - 2 * Cyr / v_x);
    alpha56 = (1/m) * (-m * v_x - 2 * Cyf * cos(delta) * Lf / v_x + 2 * Cyr * Lr / v_x);
    beta51 = (1/m) * (2 * Cyf * (cos(delta) - delta * sin(delta)) + 2 * Cxf * Sf * cos(delta));
    beta52 = (1/m) * (2 * Cxf * sin(delta));

    alpha64 = (1/Iz) * (2 * Cyf * Lf * cos(delta) * (v_y+ Lf * omega) / v_x^2 - 2 * Cyr * Lr * (v_y - Lr * omega) / v_x^2);
    alpha65 = (1/Iz) * (-2 * Cyf * Lf * cos(delta) / v_x + 2 * Cyr * Lr / v_x);
    alpha66 = (1/Iz) * (-2 * Cyf * Lf^2 * cos(delta) / v_x - 2 * Cyr * Lr^2 / v_x);
    beta61 = (1/Iz) * (2 * Cyf * Lf * (cos(delta) - delta * sin(delta)) + 2 * Cxf * Sf * Lf * cos(delta));
    beta62 = (1/Iz) * (2 * Cxf * Lf * sin(delta));

    % Recompute Ad and Bd
    A = [0 0 alpha13 alpha14 alpha15 0;
        0 0 alpha23 alpha24 alpha25 0;
        0 0 0 0 0 alpha36;
        0 0 0 alpha44 alpha45 alpha46;
        0 0 0 alpha54 alpha55 alpha56;
        0 0 0 alpha64 alpha65 alpha66];

    B = [0 0 0;
        0 0 0;
        0 0 0;
        beta41 beta42 beta43;
        beta51 beta52 0;
        beta61 beta62 0];

    Ad = eye(size(A)) + Ts * A;
    Bd = Ts * B;


    %% Controllability check (only initially)
    if k == 1 && rank(ctrb(Ad, Bd)) ~= size(Ad,1)
        error('The system is NOT controllable.');
    end

    %% LQR gain clearly computed
    K(:,:,k) = dlqr(Ad, Bd, Q, R);

    %% State deviation (error) explicitly wrapped for angles
    X_tilde(1:2, k) = X_lqr(1:2,k) - X_nom(1:2,next_index_lqr(k));
    X_tilde(3, k)   = wrapToPi(X_lqr(3,k) - X_nom(3,next_index_lqr(k)));
    X_tilde(4:6, k) = X_lqr(4:6,k) - X_nom(4:6,next_index_lqr(k));

    %% LQR input explicitly computed
    U_lqr(:, k) = U_nom(:, next_index_lqr(k)) - K(:,:,k) * X_tilde(:, k);

    %% Actuator limits
    U_lqr(1,k) = min(max(U_lqr(1,k), -delta_max), delta_max);   % steering
    U_lqr(2,k) = 0;                                             % front slip kept at zero (rear-wheel drive)
    U_lqr(3,k) = min(max(U_lqr(3,k), Sr_min), Sr_max);          % rear slip (acceleration range)
    

    % LQR Dynamic bicycle model
    dX_lqr = dynamic_bicycle_model(0, X_lqr(:,k), U_lqr(:,k), m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr);
    X_lqr(:,k+1) = X_lqr(:,k) + Ts * dX_lqr;
    
    %Limit steering within +/-Pi
    %X_lqr(3,k+1) = wrapToPi(X_lqr(3,k+1));


    % Calculate distance to the initial position to terminate iteration
    distanceToStart(k) = dis([next_lqr(1,1); next_lqr(2,1)],[X_lqr(1,k); X_lqr(2,k)]);
    if distanceToStart(k) < 5 && k > 100
        fprintf('Stopping simulation at step %d because distance is less than the threshold.\n', k);
        break;
    end
    % 
    % % Plot next waypoint and current position
    % plot(next_lqr(1,k), next_lqr(2,k), 'bo'); hold on; % Next waypoint
    % plot(X_lqr(1,k), X_lqr(2,k), 'kx');        % Current position
    % drawnow;

end

%% Performance summary
t = Tspan(1:k);
t_sat = Ts*sum(U_lqr(3,1:k) >= Sr_max - 1e-12 | U_lqr(3,1:k) <= Sr_min + 1e-12);
fprintf('\n%s\n', 'Example 7.2: TV-LQR, dynamic bicycle');
fprintf('lap time            : %.2f s\n', t(end));
fprintf('RMS CTE             : %.3f m\n', rms(e_lqr(1:k)));
fprintf('max |CTE|           : %.3f m\n', max(abs(e_lqr(1:k))));
fprintf('RMS heading error   : %.3f deg\n', rms(rad2deg(theta_e_lqr(1:k))));
fprintf('max |delta|         : %.2f deg\n', max(abs(rad2deg(U_lqr(1,1:k)))));
fprintf('time at Sr limit    : %.2f s\n', t_sat);
fprintf('RMS speed error     : %.3f m/s\n', rms(X_lqr(4,1:k) - X_nom(4,next_index_lqr(1:k))));
fprintf('speed range         : [%.1f, %.1f] m/s\n', min(X_lqr(4,1:k)), max(X_lqr(4,1:k)));
fprintf('vy range            : [%.2f, %.2f] m/s\n', min(X_lqr(5,1:k)), max(X_lqr(5,1:k)));

%% Set Global defaults--affects all figures
set(groot,'defaultFigureColor','w', ...
          'defaultAxesFontName','Helvetica', ...
          'defaultAxesFontSize',16, ...
          'defaultAxesFontWeight','bold', ...
          'defaultAxesLineWidth',1.5, ...
          'defaultLineLineWidth',3, ...
          'defaultLegendFontSize',14, ...
          'defaultLegendFontWeight','bold');
idx_k = next_index_lqr(1:k);

%% Figure 1: trajectory
figure;
plot(reference(1,:), reference(2,:), 'b', 'LineWidth', 2); hold on;
plot(X_lqr(1,1:k), X_lqr(2,1:k), 'r--', 'LineWidth', 2);
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Reference Trajectory', 'LQR Trajectory');
grid on; box on;

%% Figure 2: heading and longitudinal velocity
figure;
subplot(2,1,1); hold on; box on;
plot(t, rad2deg(X_nom(3,idx_k)));
plot(t, rad2deg(X_lqr(3,1:k)),'r--');
xlabel('Time (s)'); ylabel('Heading Angle (deg)');
legend('Nominal','LQR','Location','best'); grid on; axis("tight");
subplot(2,1,2); hold on; box on;
plot(t, velocities(idx_k), 'k:');
plot(t, X_nom(4,idx_k), 'Color', [0 0.4470 0.7410]);   % nominal in blue, as in the other plots
plot(t, X_lqr(4,1:k),'r--');
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)');
legend('Speed limit','Nominal','LQR','Location','best'); grid on; axis("tight");

%% Figure 3: lateral velocity and yaw rate
figure;
subplot(2,1,1); hold on; box on;
plot(t, X_nom(5,idx_k));
plot(t, X_lqr(5,1:k),'r--');
xlabel('Time (s)'); ylabel('Lateral Velocity (m/s)');
legend('Nominal','LQR','Location','best'); grid on; axis("tight");
subplot(2,1,2); hold on; box on;
plot(t, X_nom(6,idx_k));
plot(t, X_lqr(6,1:k),'r--');
xlabel('Time (s)'); ylabel('Yaw Rate (rad/s)');
legend('Nominal','LQR','Location','best'); grid on; axis("tight");

%% Figure 4: control inputs (steering in deg, slip ratios dimensionless)
figure;
subplot(3,1,1); hold on; box on;
plot(t, rad2deg(U_nom(1,idx_k)));
plot(t, rad2deg(U_lqr(1,1:k)),'r--');
ylabel('Steering Angle (deg)');
legend('Nominal','LQR','Location','best'); grid on; axis("tight");
subplot(3,1,2); hold on; box on;
plot(t, U_nom(2,idx_k));
plot(t, U_lqr(2,1:k),'r--');
ylabel('S_f'); grid on; axis("tight");
subplot(3,1,3); hold on; box on;
plot(t, U_nom(3,idx_k));
plot(t, U_lqr(3,1:k),'r--');
xlabel('Time (s)'); ylabel('S_r'); grid on; axis("tight");

%% Figure 5: geometric errors
figure;
subplot(2,1,1);
plot(t, e_lqr(1:k));
xlabel('Time (s)'); ylabel('CTE (m)');
grid on; box on; axis("tight");
subplot(2,1,2);
plot(t, rad2deg(theta_e_lqr(1:k)));
xlabel('Time (s)'); ylabel('Heading Error (deg)');
grid on; box on; axis("tight");

%% Save the five figures used in the book (figures/CH07)
figdir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'figures', 'CH07');
names  = {'FigEx6_2_path','FigEx6_2_states1','FigEx6_2_states2','FigEx6_2_control_input','FigEx6_2_errors'};
for i = 1:5
    exportgraphics(figure(i), fullfile(figdir, [names{i} '.jpeg']), 'Resolution', 300);
end
