%Lateral + Longitudinal control
%Minimized state error
clc; clear; close all;

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

%% Nominal inputs computed from reference (delta; Sf; Sr)
U_nom(1,1:N-1) = diff(unwrap(reference(3,:)))./Ts .* (L ./ reference(4,1:end-1)); % steering angle
U_nom(1,N) = U_nom(1,N-1);
delta_max = deg2rad(30);

Sf = 0;
U_nom(2,:) = Sf*ones(size(U_nom(1,:)));
acceleration_nominal = gradient(velocities, Ts);

Sr = (m .* acceleration_nominal) ./ (2 * Cxr);
U_nom(3,:) = Sr;
Sr_min = min(Sr); Sr_max = max(Sr); 

%% Nominal states
[x_desired, y_desired] = deal(WP(1,:), WP(2,:));
psi_desired = unwrap(headings);
vx_desired = velocities;
vy_desired = zeros(1,N); % assuming nominal lateral velocity is zero
omega_desired = gradient(psi_desired, Ts);
X_nom = [x_desired; y_desired; psi_desired; vx_desired; vy_desired; omega_desired];

%% Initial state
X_mpc = zeros(6, N);
X_mpc(:,1) = [WP(:,1); 0; 0; 0; 0];

%% MPC parameters
N_p = 15;        % Prediction horizon
N_c = 5;        % Control horizon

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
    alpha46 = (1/m) * (m * v_y - 2 * Cyf * sin(delta) * Lf / v_x);
    beta41 = (1/m) * (-2 * Cxf * Sf * sin(delta) - 2 * Cyf * (sin(delta) + delta * cos(delta)));
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

    if k ==1

        sys = ss(Ad, Bd, eye(6), 0, Ts);
        mpcobj = mpc(sys, Ts, N_p, N_c);
        mpcstate_obj = mpcstate(mpcobj);

        % MPC weights
        mpcobj.Weights.OutputVariables          = [8, 8, 5, 0.1, 0.1, 0.1];  % [x,y,yaw,vx,vy,omega]
        mpcobj.Weights.ManipulatedVariables     = [0.2, 0, 1e-3];            % steer, Sf(fixed), Sr
        mpcobj.Weights.ManipulatedVariablesRate = [3.0, 0, 1e-3];            % heavier steer rate = smoother lateral


        % Constraints
        delta_max = deg2rad(30);
        mpcobj.ManipulatedVariables(1).Min = -delta_max;
        mpcobj.ManipulatedVariables(1).Max = delta_max;
        mpcobj.ManipulatedVariables(2).Min = 0; % Constant
        mpcobj.ManipulatedVariables(2).Max = 0;
        mpcobj.ManipulatedVariables(3).Min = Sr_min;
        mpcobj.ManipulatedVariables(3).Max = Sr_max;
    end

    mpcobj.Model.Plant = ss(Ad, Bd, eye(6), zeros(6,3), Ts);

    %%%%%%%%%%%%%% Compute the control input using MPC %%%%%%%%%%%%%%%%%%%
    y = X_mpc(:,k) - X_nom(:, next_index_mpc(k)); % current state
    %y(3) = wrapToPi(y(3));
    r = zeros(size(y));     % reference output
    [U_mpc(:,k), info] = mpcmove(mpcobj, mpcstate_obj, y, r);  


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

%% Set Global defaults--affects all figures
set(groot,'defaultFigureColor','w', ...
          'defaultAxesFontName','Helvetica', ...
          'defaultAxesFontSize',16, ...
          'defaultAxesFontWeight','bold', ...
          'defaultAxesLineWidth',1.5, ...
          'defaultLineLineWidth',3, ...
          'defaultLegendFontSize',14, ...
          'defaultLegendFontWeight','bold');

%% Plot Trajectory
figure;
plot(reference(1,:), reference(2,:), 'b', 'LineWidth', 2); hold on;
plot(X_mpc(1,1:k), X_mpc(2,1:k), 'r--', 'LineWidth', 2);
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Reference Trajectory', 'mpc Trajectory');
grid on;


%% Plot States
figure;  hold on;
subplot(2,1,1);  hold on; box on;
plot(Tspan(1:k), rad2deg(X_nom(3,next_index_mpc(1:k))));
plot(Tspan(1:k), rad2deg(X_mpc(3,1:k)),'r--');
xlabel('Time (s)'); ylabel('Heading Angle (deg)');
legend('Nominal','mpc'); grid on; axis("tight");
subplot(2,1,2); hold on; box on;
plot(Tspan(1:k), X_nom(4,next_index_mpc(1:k)));
plot(Tspan(1:k), X_mpc(4,1:k),'r--');
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)');
legend('Nominal','mpc'); grid on; axis("tight");

figure;  hold on;
subplot(2,1,1);  hold on; box on;
plot(Tspan(1:k), X_nom(5,next_index_mpc(1:k)));
plot(Tspan(1:k), X_mpc(5,1:k),'r--');
xlabel('Time (s)'); ylabel('Lateral Velocity (m/s)');
legend('Nominal','mpc'); grid on; axis("tight");
subplot(2,1,2); hold on; box on;
plot(Tspan(1:k), X_nom(6,next_index_mpc(1:k)));
plot(Tspan(1:k), X_mpc(6,1:k),'r--');
xlabel('Time (s)'); ylabel('Omega (rad/s^2)');
legend('Nominal','mpc'); grid on; axis("tight");


%% Plot Control Inputs
figure; hold on;
subplot 311; 
hold on; plot(Tspan(1:k), rad2deg(U_mpc(1,1:k)));
xlabel('Time (s)'); ylabel('Steering (deg)');
grid on; box on; axis("tight");
subplot 312; 
hold on; plot(Tspan(1:k), rad2deg(U_mpc(2,1:k)));
xlabel('Time (s)'); ylabel('Sf');
grid on; box on; axis("tight");
subplot 313; 
hold on; plot(Tspan(1:k), U_mpc(3,1:k));
xlabel('Time (s)'); ylabel('Sr');
grid on; box on; axis("tight");

%% Plot Errors
figure; hold on;
subplot 211;
plot(Tspan(1:k), e_mpc(1:k));
xlabel('Time (s)'); ylabel('CTE (m)');
grid on;  box on; axis("tight");
subplot 212; 
plot(Tspan(1:k), rad2deg(theta_e_mpc(1:k)));
xlabel('Time (s)'); ylabel('Heading Error (deg)');
grid on;  box on; axis("tight");