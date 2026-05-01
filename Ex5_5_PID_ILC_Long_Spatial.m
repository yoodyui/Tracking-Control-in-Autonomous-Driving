% Example 5.5 (Spatial ILC): PID lateral + spatial ILC FF + PID speed
%  - u_pid(k): nominal lateral PID steering
%  - uL_s(i): learned steering FF indexed by waypoint i (spatial, not time)
%  - delta(k) = u_pid(k) + uL_s(i_k)
%  - Speed is controlled by a separate PID (no learning)
% Don't need an extra ZPF to stabilize the system
clear; clc; close all;

%% Parameters
L       = 1.08;                % Wheelbase [m]
Ts      = 0.05;                % Sampling time [s]
N       = 583;                 % Max steps per traversal (can stop earlier)
n_iter  = 101;                  % Number of iterations (laps)
del_max = deg2rad(30);         % Steering saturation [rad]

%% Reference (same path each iteration)
n_wp = 100;
[WP, velocities, headings] = waypoints_with_velocity_and_heading(n_wp);
reference = [WP; unwrap(headings); velocities];

%% Lateral PID (for steering u_pid)
Kp = 0.10;  Ki = 0.00;  Kd = 0.00;  Kt = 1.00;   % tune as needed

%% Longitudinal PID (for speed -> acceleration u(1,k))
KvP = 0.8;  KvI = 0.2;  KvD = 0.05;              % tune as needed
a_min = -3.0;   % [m/s^2] decel limit
a_max =  2.0;   % [m/s^2] accel limit

%% Spatial ILC (lateral FF only): uL_s^{j+1}(i) = uL_s^j(i) + L_ilc * e_avg(i)
L_ilc = 0.01;

%% Storage across iterations
x_all    = cell(1, n_iter);            % states per iteration
e_all    = cell(1, n_iter);            % CTE per iteration (time-indexed for plotting)
rmse_cte = zeros(1, n_iter);

%% Initialize learned lateral FF (spatial, one value per waypoint)
uL_s = zeros(1, n_wp);                 % spatial FF steering per waypoint

for j = 1:n_iter

    % State: [x; y; psi; vx]
    x = zeros(4, N+1);
    x(:,1) = [WP(1:2,1); 0; 0];        % [x; y; psi; vx]

    % Inputs per time step:
    %   u(1,k) = longitudinal acceleration (from speed PID)
    %   u(2,k) = steering = u_pid (lateral) + uL_s(i_k) (spatial FF)
    u = zeros(2, N);

    % Lateral PID integrator and book-keeping
    int_cte = 0;

    % Speed PID integrator and book-keeping
    iv = 0;                  % integral of speed error
    ev_prev = 0;             % previous speed error

    % Logs (time-indexed for this lap)
    e_cte    = zeros(1, N);  % CTE
    theta_e  = zeros(1, N);  % heading error
    ev       = zeros(1, N);  % speed error
    idx_used = zeros(1, N);  % waypoint indices used each step

    % Spatial accumulators (for ILC update)
    e_sum = zeros(1, n_wp);
    e_cnt = zeros(1, n_wp);

    for k = 1:N
        % Nearest waypoint and segment for CTE computation
        [next_k, pass_k, i_k, ~] = nextWP(WP, x(1:2,k));  
        idx_used(k) = i_k;

        % Desired states from nearest waypoint
        Xr   = reference(1, i_k);
        Yr   = reference(2, i_k);
        psir = reference(3, i_k);
        vref = reference(4, i_k);

        % Lateral errors
        e_cte(k)   = cte(next_k, pass_k, x(1:2,k));
        theta_e(k) = atan2( sin(psir - x(3,k)), cos(psir - x(3,k)) );

        % ---------- Longitudinal PID (speed control) ----------
        ev(k)   = vref - x(4,k);                         % speed error
        dev     = (k==1)*0 + (k>1)*(ev(k) - ev_prev)/Ts; % derivative
        iv_try  = iv + ev(k)*Ts;                         % tentative I
        a_unsat = KvP*ev(k) + KvI*iv_try + KvD*dev;      % accel command (unsat)
        a_sat   = min(max(a_unsat, a_min), a_max);       % saturate
        if abs(a_unsat - a_sat) < 1e-9
            iv = iv_try;                                 % anti-windup (simple)
        end
        u(1,k)  = a_sat;                                 % apply accel
        ev_prev = ev(k);

        % ---------- Lateral PID (nominal steering) ----------
        if k == 1
            u_pid = Kp*e_cte(k) + Ki*int_cte + Kt*theta_e(k);
        else
            de_cte  = (e_cte(k) - e_cte(k-1))/Ts;
            int_cte = int_cte + e_cte(k)*Ts;
            u_pid   = Kp*e_cte(k) + Ki*int_cte + Kd*de_cte + Kt*theta_e(k);
        end

        % ---------- Spatial FF: pick value at waypoint i_k ----------
        delta = u_pid + uL_s(i_k);
        delta = max(min(delta, del_max), -del_max);
        u(2,k) = delta;

        % ---------- Kinematic bicycle update ----------
        x(:,k+1) = x(:,k) + Ts * [ ...
            x(4,k)*cos(x(3,k));
            x(4,k)*sin(x(3,k));
            (x(4,k)/L)*tan(delta);
            u(1,k) ];

        % Accumulate spatial error at this waypoint (for ILC update)
        e_sum(i_k) = e_sum(i_k) + e_cte(k);
        e_cnt(i_k) = e_cnt(i_k) + 1;

        % Early stop if back near start
        if k > 100
            d0 = hypot(WP(1,1) - x(1,k), WP(2,1) - x(2,k));
            if d0 < 5
                fprintf('Iter %d: stopping at step %d (near start)\n', j, k);
                % Truncate logs to actual length k
                e_cte    = e_cte(1:k);
                theta_e  = theta_e(1:k);
                ev       = ev(1:k);
                idx_used = idx_used(1:k);
                u        = u(:,1:k);
                x        = x(:,1:k+1);
                break;
            end
        end
    end

    % Save iteration data (time-indexed traces)
    x_all{j} = x;
    e_all{j} = e_cte;
    rmse_cte(j) = sqrt(mean(e_cte.^2));

    % ---------- Spatial ILC update ----------
    e_avg = e_sum ./ max(1, e_cnt);     % average CTE per waypoint visited
    e_avg(~isfinite(e_avg)) = 0;        % guard (in case of no visits)
    uL_s = uL_s + L_ilc * e_avg;        % update spatial FF

    % Optional: mild smoothing around the loop
    uL_s = smoothdata(uL_s, 'movmean', 5);

    % Wrap-around continuity (optional, for a closed loop)
    % enforce end≈start continuity by averaging ends
    uL_s(1) = 0.5*(uL_s(1) + uL_s(end));
    uL_s(end) = uL_s(1);
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
title('Spatial ILC: PID lateral + spatial FF with PID speed');
lgd = [{'Reference'}, arrayfun(@(jj) sprintf('Iter %d', jj), 1:n_iter, 'UniformOutput', false)];
legend(lgd{:}, 'Location', 'best');

%% 3D plot of CTE across iterations (time-indexed view)
figure; hold on; grid on; box on;
for j = 1:n_iter
    e = e_all{j};
    t = (0:numel(e)-1)*Ts;
    plot3(t, (j-1)*ones(size(t)), e, 'LineWidth', 1.3);
end
xlabel('Time t (s)'); ylabel('Iteration j'); zlabel('CTE e_j(k) (m)');
title('CTE Trajectories Across Iterations'); view(135,25);

%% RMSE vs iteration
figure;
plot(0:n_iter-1, rmse_cte, '-o', 'LineWidth', 1.8, 'MarkerSize', 5);
grid on; box on;
xlabel('Iteration j'); ylabel('RMSE(e_j) (m)');
title('CTE RMSE per Iteration (Spatial ILC)');

%% (Optional) Spatial FF profile learned
figure; 
plot(1:n_wp, rad2deg(uL_s), 'LineWidth', 1.8);
grid on; box on; xlabel('Waypoint index i'); ylabel('u_L^s(i) (deg)');
title('Learned Spatial Steering Feedforward');
