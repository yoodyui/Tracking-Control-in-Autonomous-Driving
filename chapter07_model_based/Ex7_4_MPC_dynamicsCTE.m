% Example 7.4: Path-Indexed MPC with Output/Path-Error Tracking
%Lateral + Longitudinal control
%Minimize the output [x,y,cte,psi,vx,vy,omega]
clc; clear; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

%% Vehicle Parameters
m = 1500;       % Vehicle mass (kg)
Lf = 1.2;       % Distance from CG to front axle (m)
Lr = 1.6;       % Distance from CG to rear axle (m)
L = Lf + Lr;
Iz = m * ((0.5 * (Lf + Lr))^2);  % Yaw moment of inertia (kg*m^2)
Cdrag = 0.32;      % Drag coefficient
A = 2.2;        % Frontal area (m^2)
rho = 1.225;    % Air density (kg/m^3)
CxA = 0.5 * Cdrag * A * rho;
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

%% Initial state
X_mpc = zeros(6, N);
X_mpc(:,1) = [WP(:,1); psi_desired(1); vx_desired(1); 0; 0];   % start on the path at the nominal speed

%% MPC parameters
N_p = 15;        % Prediction horizon
N_c = 5;        % Control horizon
delta_max = deg2rad(30);

%% Simulation Loop
for k = 1:N-1
   
        %% Nearest forward waypoint and geometric errors
    % if k == 1
    %     [next_mpc(:,k),pass_mpc(:,k),next_index_mpc(k),pass_index_mpc(k)] = nextWP(WP,[X_mpc(1,k); X_mpc(2,k)]);
    % else
    %     [next_mpc(:,k),pass_mpc(:,k),next_index_mpc(k),pass_index_mpc(k)] = nextForwardWP(WP,[X_mpc(1,k); X_mpc(2,k)], next_index_mpc(k-1), 30); % window=30 (tune 20–40)
    % end
    [next_mpc(:,k),pass_mpc(:,k),next_index_mpc(k),pass_index_mpc(k)] = nextWP(WP,[X_mpc(1,k); X_mpc(2,k)]);


    % Compute lateral error from mpc
    e_mpc(k) = cte(next_mpc(:,k),pass_mpc(:,k),[X_mpc(1,k); X_mpc(2,k)]);
    %Calculate heading error from mpc
    theta_e_mpc(k) = mod(headings(next_index_mpc(k)) - X_mpc(3,k) + pi, 2*pi) - pi;

    %% Compute Longitudinal error from mpc along heading
    wp_error_vector = WP(:,next_index_mpc(k)) - [X_mpc(1,k); X_mpc(2,k)];
    long_error(k) = dot(wp_error_vector, [cos(X_mpc(3,k)); sin(X_mpc(3,k))]);

    % Update nominal states
    x = X_nom(1,next_index_mpc(k));
    y = X_nom(2,next_index_mpc(k));
    psi = X_nom(3,next_index_mpc(k));
    v_x = X_nom(4,next_index_mpc(k));
    v_y = X_nom(5,next_index_mpc(k));
    omega = X_nom(6,next_index_mpc(k));

    %control inputs
    delta = U_nom(1,next_index_mpc(k));

   
    %% mpc design
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

    % --- Output mapping: y = Ck * x  ---
    % cte = -sin(psi)*(x-xr) + cos(psi)*(y-yr)
    % Implement as C*x and put the reference shift into r(1) = cte_ref = 0
    Cd = [1 0 0 0 0 0;    % x
        0 1 0 0 0 0;    % y
        -sin(psi) cos(psi) 0 0 0 0;  % cte (optional)
        0 0 1 0 0 0;    % yaw
        0 0 0 1 0 0;   % vx
        0 0 0 0 1 0;    %vy
        0 0 0 0 0 1];   %omega

    Dd = zeros(length(Cd),3);


    if k ==1

        sys = ss(Ad, Bd, Cd, Dd, Ts);
        mpcobj = mpc(sys, Ts, N_p, N_c);
        mpcstate_obj = mpcstate(mpcobj);

        % MPC weights
        mpcobj.Weights.OutputVariables = [1e-4, 1e-4, 8, 3, 0.5, 1e-4, 1e-4];  %x, y, cte, yaw, vx, vy, omega
        mpcobj.Weights.ManipulatedVariables     = [0.2, 0, 1e-4];            % steer, Sf(fixed), Sr
        mpcobj.Weights.ManipulatedVariablesRate = [3.0, 0, 1e-3];            % heavier steer rate = smoother lateral


        % Constraints. The MPC optimizes the input deviation u_tilde = u - u_nom,
        % so the widest possible deviation range is declared here and the exact
        % bounds u_min - u_nom(i_k) <= u_tilde <= u_max - u_nom(i_k) are passed
        % at run time through mpcmoveopt.
        mpcobj.ManipulatedVariables(1).Min = -2*delta_max;
        mpcobj.ManipulatedVariables(1).Max =  2*delta_max;
        mpcobj.ManipulatedVariables(2).Min = 0; % front slip fixed at zero
        mpcobj.ManipulatedVariables(2).Max = 0;
        mpcobj.ManipulatedVariables(3).Min = Sr_min - Sr_max;
        mpcobj.ManipulatedVariables(3).Max = Sr_max - Sr_min;
        mpcopt = mpcmoveopt;
    end

    mpcobj.Model.Plant = ss(Ad, Bd, Cd, Dd, Ts);

    %%%%%%%%%%%%%% Compute the control input using MPC %%%%%%%%%%%%%%%%%%%
    % Outputs are formed from the state deviation, y = Cd*(x - x_nom(i_k)),
    % so that the third output is the cross-track error and every reference is zero.
    u_nom_k = U_nom(:, next_index_mpc(k));              % nominal input at the nearest waypoint
    x_til = X_mpc(:,k) - X_nom(:, next_index_mpc(k));   % state deviation
    x_til(3) = wrapToPi(x_til(3));
    y = Cd * x_til;                                     % [x; y; cte; yaw; vx; vy; omega] deviations
    r = zeros(size(y));
    mpcopt.MVMin = ([-delta_max; 0; Sr_min] - u_nom_k).';   % bounds on u_tilde
    mpcopt.MVMax = ([ delta_max; 0; Sr_max] - u_nom_k).';
    [u_tilde, info] = mpcmove(mpcobj, mpcstate_obj, y, r, [], mpcopt);
    U_mpc(:,k) = u_nom_k + u_tilde;                     % applied input u = u_nom + u_tilde


    % MPC Dynamic bicycle model
    dX_mpc = dynamic_bicycle_model(0, X_mpc(:,k), U_mpc(:,k), m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr);
    X_mpc(:,k+1) = X_mpc(:,k) + Ts * dX_mpc;
    
    %Limit steering within +/-Pi
    %X_mpc(3,k+1) = wrapToPi(X_mpc(3,k+1));


    % Calculate distance to the initial position to terminate iteration
    distanceToStart(k) = dis([next_mpc(1,1); next_mpc(2,1)],[X_mpc(1,k); X_mpc(2,k)]);
    if distanceToStart(k) < 5 && k > 100
        fprintf('Stopping simulation at step %d because distance is less than the threshold.\n', k);
        break;
    end

    % % Plot next waypoint and current position
    % plot(next_mpc(1,k), next_mpc(2,k), 'bo'); hold on; % Next waypoint
    % plot(X_mpc(1,k), X_mpc(2,k), 'kx');        % Current position
    % drawnow;

end

%% Performance summary
t = Tspan(1:k);
t_sat = Ts*sum(U_mpc(3,1:k) >= Sr_max - 1e-12 | U_mpc(3,1:k) <= Sr_min + 1e-12);
fprintf('\n%s\n', 'Example 7.4: MPC (output/path-error), dynamic bicycle');
fprintf('lap time            : %.2f s\n', t(end));
fprintf('RMS CTE             : %.3f m\n', rms(e_mpc(1:k)));
fprintf('max |CTE|           : %.3f m\n', max(abs(e_mpc(1:k))));
fprintf('RMS heading error   : %.3f deg\n', rms(rad2deg(theta_e_mpc(1:k))));
fprintf('max |delta|         : %.2f deg\n', max(abs(rad2deg(U_mpc(1,1:k)))));
fprintf('time at Sr limit    : %.2f s\n', t_sat);
fprintf('RMS speed error     : %.3f m/s\n', rms(X_mpc(4,1:k) - X_nom(4,next_index_mpc(1:k))));
fprintf('speed range         : [%.1f, %.1f] m/s\n', min(X_mpc(4,1:k)), max(X_mpc(4,1:k)));
fprintf('vy range            : [%.2f, %.2f] m/s\n', min(X_mpc(5,1:k)), max(X_mpc(5,1:k)));

%% Set Global defaults--affects all figures
set(groot,'defaultFigureColor','w', ...
          'defaultAxesFontName','Helvetica', ...
          'defaultAxesFontSize',16, ...
          'defaultAxesFontWeight','bold', ...
          'defaultAxesLineWidth',1.5, ...
          'defaultLineLineWidth',3, ...
          'defaultLegendFontSize',14, ...
          'defaultLegendFontWeight','bold');
idx_k = next_index_mpc(1:k);

%% Figure 1: trajectory
figure;
plot(reference(1,:), reference(2,:), 'b', 'LineWidth', 2); hold on;
plot(X_mpc(1,1:k), X_mpc(2,1:k), 'r--', 'LineWidth', 2);
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Reference Trajectory', 'MPC Trajectory');
grid on; box on;

%% Figure 2: heading and longitudinal velocity
figure;
subplot(2,1,1); hold on; box on;
plot(t, rad2deg(X_nom(3,idx_k)));
plot(t, rad2deg(X_mpc(3,1:k)),'r--');
xlabel('Time (s)'); ylabel('Heading Angle (deg)');
legend('Nominal','MPC','Location','best'); grid on; axis("tight");
subplot(2,1,2); hold on; box on;
plot(t, velocities(idx_k), 'k:');
plot(t, X_nom(4,idx_k), 'Color', [0 0.4470 0.7410]);   % nominal in blue, as in the other plots
plot(t, X_mpc(4,1:k),'r--');
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)');
legend('Speed limit','Nominal','MPC','Location','best'); grid on; axis("tight");

%% Figure 3: lateral velocity and yaw rate
figure;
subplot(2,1,1); hold on; box on;
plot(t, X_nom(5,idx_k));
plot(t, X_mpc(5,1:k),'r--');
xlabel('Time (s)'); ylabel('Lateral Velocity (m/s)');
legend('Nominal','MPC','Location','best'); grid on; axis("tight");
subplot(2,1,2); hold on; box on;
plot(t, X_nom(6,idx_k));
plot(t, X_mpc(6,1:k),'r--');
xlabel('Time (s)'); ylabel('Yaw Rate (rad/s)');
legend('Nominal','MPC','Location','best'); grid on; axis("tight");

%% Figure 4: control inputs (steering in deg, slip ratios dimensionless)
figure;
subplot(3,1,1); hold on; box on;
plot(t, rad2deg(U_nom(1,idx_k)));
plot(t, rad2deg(U_mpc(1,1:k)),'r--');
ylabel('Steering Angle (deg)');
legend('Nominal','MPC','Location','best'); grid on; axis("tight");
subplot(3,1,2); hold on; box on;
plot(t, U_nom(2,idx_k));
plot(t, U_mpc(2,1:k),'r--');
ylabel('S_f'); grid on; axis("tight");
subplot(3,1,3); hold on; box on;
plot(t, U_nom(3,idx_k));
plot(t, U_mpc(3,1:k),'r--');
xlabel('Time (s)'); ylabel('S_r'); grid on; axis("tight");

%% Figure 5: geometric errors
figure;
subplot(2,1,1);
plot(t, e_mpc(1:k));
xlabel('Time (s)'); ylabel('CTE (m)');
grid on; box on; axis("tight");
subplot(2,1,2);
plot(t, rad2deg(theta_e_mpc(1:k)));
xlabel('Time (s)'); ylabel('Heading Error (deg)');
grid on; box on; axis("tight");

%% Save the five figures used in the book (figures/CH07)
figdir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'figures', 'CH07');
names  = {'FigEx6_4_path','FigEx6_4_states1','FigEx6_4_states2','FigEx6_4_control_input','FigEx6_4_errors'};
for i = 1:5
    exportgraphics(figure(i), fullfile(figdir, [names{i} '.jpeg']), 'Resolution', 300);
end
