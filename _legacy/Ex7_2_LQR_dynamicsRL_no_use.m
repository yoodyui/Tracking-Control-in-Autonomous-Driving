% Example 7.2 auxiliary script: RL-tuned LQR variant (not used in text)
%Picking the closest WP in the loop
%Work very well for all velocity!!
%Add RL to adjust Q and R

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

%% Nominal inputs (compute from reference)
U_nom(1,1:N-1) = diff(unwrap(reference(3,:)))./Ts .* (L ./ reference(4,1:end-1)); % steering angle
U_nom(1,N) = U_nom(1,N-1);

Sf = 0;
U_nom(2,:) = Sf*ones(size(U_nom(1,:)));

acceleration_nominal = gradient(velocities, Ts);
Sr = (m .* acceleration_nominal) ./ (2 * Cxr);
U_nom(3,:) = Sr;

%% Find the nominal path
X0_nominal = [WP(:,1); headings(1); velocities(1); 0; 0];

[t_nl, X_nl] = ode45(@(t, X) dynamic_bicycle_model(t, X, ...
    [interp1(Tspan, U_nom(1,:), t, 'previous', 'extrap'); 
     Sf;  % Sf is constant
     interp1(Tspan, Sr, t, 'previous', 'extrap')], ... % explicitly interpolate Sr
    m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr), Tspan, X0_nominal);

X_nl = X_nl';


%% Nominal states
[x_desired, y_desired] = deal(WP(1,:), WP(2,:));
psi_desired = headings;
vx_desired = velocities;
vy_desired = zeros(1,N); % assuming nominal lateral velocity is zero
omega_desired = gradient(unwrap(psi_desired), Ts);
X_nom = [x_desired; y_desired; psi_desired; vx_desired; vy_desired; omega_desired];

%% Initialize states
X_lqr = zeros(6, N);
X_lqr(:,1) = [WP(:,1); headings(1); velocities(1); 0; 0];

%% LQR Weight Matrices
Q = diag([1, 1, 1, 1, 0.01, 1]);
R = diag([1, 1, 1]);

%% Simple RL-based Adaptive Q and R updates
alpha = 3; % learning rate for Q
beta = 0.1;  % learning rate for R

% Desired errors (targets)
target_error = 0; %desired values are zero
delta_target = deg2rad(1); % very small target steering

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

    %% Explicit saturation per input clearly
    % U_lqr(1,k) = min(max(U_lqr(1,k), deg2rad(-30)), deg2rad(30)); % steering clearly saturated
    % U_lqr(2,k) = 0; % consistent with zero front slip assumption
    % U_lqr(3,k) = min(max(U_lqr(3,k), 0.1*Sr(k)), 0.1*Sr(k)); % rear slip clearly limited
    


    % LQR Dynamic bicycle model
    dX_lqr = dynamic_bicycle_model(0, X_lqr(:,k), U_lqr(:,k), m, Iz, Lf, Lr, CxA, Cxf, Cxr, Cyf, Cyr);
    X_lqr(:,k+1) = X_lqr(:,k) + Ts * dX_lqr;
    
    %Limit steering within +/-Pi
    X_lqr(3,k+1) = wrapToPi(X_lqr(3,k+1));

    % % Velocity and yaw-rate errors explicitly
    % vx_error(k) = X_lqr(4,k) - v_x;
    % vy_error(k) = X_lqr(5,k) - v_y;
    % omega_error(k) = X_lqr(6,k) - omega;

    %Adaptive update for Q (based on tracking errors)
    Q(1,1) = Q(1,1) + alpha*(long_error(k)^2 - target_error^2);
    Q(2,2) = Q(2,2) + alpha*(e_lqr(k)^2 - target_error^2);
    Q(3,3) = Q(3,3) + alpha*(theta_e_lqr(k)^2 - target_error^2);

    % %Adaptive updates based on velocity errors (add these terms poorer
    % performance)
    % Q(4,4) = Q(4,4) + alpha*(vx_error(k)^2 - target_error^2);
    % Q(5,5) = Q(5,5) + alpha*(vy_error(k)^2 - target_error^2);
    % Q(6,6) = Q(6,6) + alpha*(omega_error(k)^2 - target_error^2);

    q(:,:,k) = Q; 

    % Adaptive updates based on R
    R(1,1) = R(1,1) + beta*(delta^2 - delta_target^2);

    r(:,:,k) = R; 

    % Calculate distance to the initial position to terminate iteration
    distanceToStart(k) = dis([next_lqr(1,1); next_lqr(2,1)],[X_lqr(1,k); X_lqr(2,k)]);
    if distanceToStart(k) < 7 && k > 100
        fprintf('Stopping simulation at step %d because distance is less than the threshold.\n', k);
        break;
    end

    % Plot next waypoint and current position
    plot(next_lqr(1,k), next_lqr(2,k), 'bo'); hold on; % Next waypoint
    plot(X_lqr(1,k), X_lqr(2,k), 'kx');        % Current position
    drawnow;

end


%% Plot Trajectory
figure;
plot(reference(1,:), reference(2,:), 'b:', 'LineWidth', 2); hold on;
%plot(X_nl(1,:), X_nl(2,:), 'r-.', 'LineWidth', 2);
plot(X_lqr(1,1:k), X_lqr(2,1:k), 'k--', 'LineWidth', 2);
xlabel('X Position (m)');
ylabel('Y Position (m)');
%legend('Reference Trajectory', 'Nominal Trajectory', 'LQR Trajectory');
legend('Reference Trajectory', 'LQR Trajectory');
grid on;
axis equal;


%% Plot States
figure;  hold on;
subplot(2,1,1);  hold on;
plot(Tspan(1:k), X_nom(3,next_index_lqr(1:k)));
plot(Tspan(1:k), X_lqr(3,1:k),'r--');
xlabel('Time (s)'); ylabel('Yaw Angle (rad)');
legend('Nominal','LQR'); grid on;
subplot(2,1,2); hold on;
plot(Tspan(1:k), X_nom(4,next_index_lqr(1:k)));
plot(Tspan(1:k), X_lqr(4,1:k),'r--');
xlabel('Time (s)'); ylabel('Longitudinal Velocity (m/s)');
legend('Nominal','LQR'); grid on;

figure;  hold on;
subplot(2,1,1);  hold on;
plot(Tspan(1:k), X_nom(5,next_index_lqr(1:k)));
plot(Tspan(1:k), X_lqr(5,1:k),'r--');
xlabel('Time (s)'); ylabel('Lateral Velocity (m/s)');
legend('Nominal','LQR'); grid on;
subplot(2,1,2); hold on;
plot(Tspan(1:k), X_nom(6,next_index_lqr(1:k)));
plot(Tspan(1:k), X_lqr(6,1:k),'r--');
xlabel('Time (s)'); ylabel('Omega (rad/s^2)');
legend('Nominal','LQR'); grid on;


%% Plot Control Inputs
figure; hold on;
subplot 311; plot(Tspan(1:k), rad2deg(U_nom(1,next_index_lqr(1:k))));
hold on; plot(Tspan(1:k), rad2deg(U_lqr(1,1:k)),'r--');
xlabel('Time (s)'); ylabel('Steering Angle (deg)');
legend('Nominal','LQR'); grid on; box on;
subplot 312; plot(Tspan(1:k), rad2deg(U_nom(2,next_index_lqr(1:k))));
hold on; plot(Tspan(1:k), rad2deg(U_lqr(2,1:k)),'r--');
xlabel('Time (s)'); ylabel('Sf');
legend('Nominal','LQR'); grid on; box on;
subplot 313; plot(Tspan(1:k), rad2deg(U_nom(3,next_index_lqr(1:k))));
hold on; plot(Tspan(1:k), rad2deg(U_lqr(3,1:k)),'r--');
xlabel('Time (s)'); ylabel('Sr');
legend('Nominal','LQR'); grid on; box on;

%% Plot Errors
figure; hold on;
subplot 311;
plot(Tspan(1:k), e_lqr(1:k));
xlabel('Time (s)'); ylabel('CTE (m)');
grid on;  box on;
subplot 312; 
plot(Tspan(1:k), theta_e_lqr(1:k));
xlabel('Time (s)'); ylabel('Heading Error (m)');
grid on;  box on;
subplot 313; plot(Tspan(1:k), long_error(1:k));
xlabel('Time (s)'); ylabel('Longitudinal Error (m)'); grid on; box on;

% % After simulation loop:
% figure;
% subplot(2,1,1);
% plot(diag(Q), 'o-');
% title('Adapted Q diagonal weights');
% xlabel('State index'); ylabel('Value'); grid on;
% 
% subplot(2,1,2);
% plot(diag(R), 's-');
% title('Adapted R diagonal weights');
% xlabel('Control input index'); ylabel('Value'); grid on;
q1(1:k) = (q(1,1,1:k));
q2(1:k) = (q(2,2,1:k));
q3(1:k) = (q(3,3,1:k));
figure;plot(q1); hold on; plot(q2,'r--'); plot(q3,'k:');
