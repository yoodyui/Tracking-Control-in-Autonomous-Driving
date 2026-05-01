% Example 5.5: PID + ILC (lateral-only learning) + PID speed control
% u_j(k): nominal PID steering, uL_j(k): learned lateral FF,
% delta_j(k) = u_j(k) + uL_j(k), e_j(k): cross-track error
% Need a ZPF to stabilize the system
clear; clc; close all;

%% Parameters
L       = 1.08;                % Wheelbase [m]
Ts      = 0.05;                % Sampling time [s]
N       = 583;                 % Steps per traversal
n_iter  = 101;                  % Number of iterations (laps)
del_max = deg2rad(30);         % Steering saturation [rad]

%% Reference (same path each iteration)
n_wp = 100;
[WP, velocities, headings] = waypoints_with_velocity_and_heading(n_wp);
reference = [WP; unwrap(headings); velocities];

%% Lateral PID (for steering u_j)
Kp = 0.10;  Ki = 0.00;  Kd = 0.00;  Kt = 1.00;   % tune as needed

%% Longitudinal PID (for speed -> acceleration u(1,k))
KvP = 0.8;  KvI = 0.2;  KvD = 0.05;            % tune as needed
a_min = -3.0;   % [m/s^2] decel limit
a_max =  2.0;   % [m/s^2] accel limit

%% ILC update (lateral FF only): uL_{j+1}(k) = uL_j(k) + L_ilc * e_j(k)
L_ilc = 0.01;

%% Storage across iterations
x_all = cell(1, n_iter);        % states per iteration
e_all = cell(1, n_iter);        % CTE per iteration
rmse_cte = zeros(1, n_iter);

%% Initialize learned lateral FF for iteration 1
uL = zeros(1, N);               % u^L_1(k) = 0

for j = 1:n_iter
    % State: [x; y; psi; vx]
    x = zeros(4, N+1);
    x(:,1) = [WP(1:2,1); 0; 0];  % [x; y; psi (heading); vx]

    % Inputs per time step:
    %   u(1,k) = longitudinal acceleration (from speed PID)
    %   u(2,k) = steering = u_pid (lateral) + uL (FF)
    u = zeros(2, N);

    % Lateral PID integrator and book-keeping
    int_cte = 0;

    % Speed PID integrator and book-keeping
    iv = 0;            % integral of speed error
    ev_prev = 0;       % previous speed error

    % Logs
    e       = zeros(1, N);     % CTE
    theta_e = zeros(1, N);     % heading error
    ev      = zeros(1, N);     % speed error

    for k = 1:N
        % Nearest waypoint and segment for CTE computation
        [next(:,k), pass(:,k), next_index(k), pass_index(k)] = nextWP(WP, [x(1,k); x(2,k)]);

        % Desired states from nearest waypoint
        x_desired(1,k) = reference(1, next_index(k));   % X_ref
        x_desired(2,k) = reference(2, next_index(k));   % Y_ref
        x_desired(3,k) = reference(3, next_index(k));   % psi_ref
        x_desired(4,k) = reference(4, next_index(k));   % v_ref

        % Lateral errors
        e_cte(k) = cte(next(:,k), pass(:,k), [x(1,k); x(2,k)]);
        theta_e(k) = atan2( sin(x_desired(3,k) - x(3,k)), cos(x_desired(3,k) - x(3,k)) );

        % ---------- Longitudinal PID (speed control) ----------
        v_ref   = x_desired(4,k);
        ev(k)   = v_ref - x(4,k);                     % speed error
        dev     = (k==1)*0 + (k>1)*(ev(k)-ev_prev)/Ts; % derivative
        iv_try  = iv + ev(k)*Ts;                      % tentative I
        a_unsat = KvP*ev(k) + KvI*iv_try + KvD*dev;   % accel command (unsat)
        a_sat   = min(max(a_unsat, a_min), a_max);    % saturate
        if abs(a_unsat - a_sat) < 1e-9
            iv = iv_try;                              % anti-windup (simple)
        end
        u(1,k)  = a_sat;                              % apply accel
        ev_prev = ev(k);

        % ---------- Lateral PID (nominal steering) ----------
        if k == 1
            u_pid = Kp*e_cte(k) + Ki*int_cte + Kt*theta_e(k);
        else
            de_cte  = (e_cte(k) - e_cte(k-1))/Ts;
            int_cte = int_cte + e_cte(k)*Ts;
            u_pid   = Kp*e_cte(k) + Ki*int_cte + Kd*de_cte + Kt*theta_e(k);
        end

        % Applied steering = nominal PID + learned FF (lateral only)
        delta = u_pid + uL(k);
        delta = max(min(delta, del_max), -del_max);
        u(2,k) = delta;

        % ---------- Kinematic bicycle update ----------
        x(:,k+1) = x(:,k) + Ts * [ ...
            x(4,k)*cos(x(3,k));
            x(4,k)*sin(x(3,k));
            (x(4,k)/L)*tan(delta);
            u(1,k) ];

        % Early stop if back near start
        distanceToStart = hypot(WP(1,1)-x(1,k), WP(2,1)-x(2,k));
        if distanceToStart < 5 && k > 100
            fprintf('Iter %d: stopping at step %d (near start)\n', j, k);
            break;
        end
    end

    % Save iteration data
    x_all{j} = x;
    e_all{j} = e_cte;
    rmse_cte(j) = sqrt(mean(e_cte.^2));   % RMSE for iteration j

    % ---------- ILC update (lateral FF only) ----------
    % Learn from the CTE sequence (same horizon N)
    uL = uL + L_ilc * e_cte;
end

%% Plots: path tracking across iterations
figure; hold on; grid on; box on;
plot(WP(1,:), WP(2,:), 'k--', 'LineWidth', 2);
C = lines(n_iter);
for j = 1:n_iter
    xx = x_all{j};
    plot(xx(1,:), xx(2,:), 'Color', C(j,:), 'LineWidth', 1.5);
end
xlabel('X (m)'); ylabel('Y (m)');
title('PID + ILC (lateral FF) with PID Speed Control');
lgd = [{'Reference'}, arrayfun(@(jj) sprintf('Iter %d', jj), 1:n_iter, 'UniformOutput', false)];
legend(lgd{:}, 'Location', 'best');

%% 3D plot of CTE across iterations
figure; hold on; grid on; box on;
for j = 1:n_iter
    e = e_all{j};
    t = (0:numel(e)-1)*Ts;
    plot3(t, (j-1)*ones(size(t)), e, 'LineWidth', 1.5);
end
xlabel('Time t (s)'); ylabel('Iteration j'); zlabel('CTE e_j(k) (m)');
title('CTE Trajectories Across Iterations'); view(135,25);

%% RMSE vs iteration
figure; 
plot(0:n_iter-1, rmse_cte, '-o', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on; box on; xlabel('Iteration j'); ylabel('RMSE(e_j) (m)');
title('CTE RMSE per Iteration (PID+ILC)');

%% (Optional) Show speed tracking on last iteration
j = n_iter; xx = x_all{j};
% Reconstruct the nearest indices for plotting (simple pass)
idx_last = zeros(1, N);
for k = 1:N
    [~,~,idx_last(k),~] = nextWP(WP, xx(1:2,k));
end
v_ref_last = reference(4, idx_last);
t = (0:N-1)*Ts;
figure; hold on; grid on; box on;
plot(t, v_ref_last, 'b--', 'LineWidth', 1.8);
plot(t, xx(4,1:N), 'r:', 'LineWidth', 1.8);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('v_{ref}','v_x','Location','best');
title(sprintf('Speed Tracking (Iteration %d)', j));
