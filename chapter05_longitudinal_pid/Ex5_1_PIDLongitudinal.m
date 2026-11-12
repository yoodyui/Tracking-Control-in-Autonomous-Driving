% Example 5.1: both cases from the revised textbook listing.
% Example: PID speed control on a straight path (delta_k = 0)
clear; clc; close all;
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'setup_paths.m'));

% Vehicle and simulation parameters
L  = 1.08;                  % Wheelbase [m]
Ts = 0.1;                   % Sampling time [s]
N  = 300;                   % Simulation steps
t  = (0:N-1)*Ts;

% PID gains
gains.Kp  = 0.9;
gains.Ki  = 0.25;
gains.Kd  = 0.05;
gains.Kff = 1.0;

% Trapezoidal reference: accelerate -> cruise -> decelerate
vCruise = 6.0;              % [m/s]
alpha   = 0.6;              % [m/s^2]
tRise   = vCruise/alpha;
tCruise = 10;               % [s]
k1 = round(tRise/Ts);
k2 = k1 + round(tCruise/Ts);

vref = zeros(1,N);
for k = 1:N
    if k <= k1
        vref(k) = min(alpha*t(k), vCruise);
    elseif k <= k2
        vref(k) = vCruise;
    else
        vref(k) = max(vCruise - alpha*(t(k)-t(k2)), 0);
    end
end

% Case A: no acceleration limit.
caseA = simulateSpeedPID(vref, Ts, gains, L, [-Inf, Inf]);

% Case B: tight acceleration limit.
caseB = simulateSpeedPID(vref, Ts, gains, L, [-0.5, 0.5]);

% Plot the two cases.
figure; plot(t,vref,'k--',t,caseA.v,'b','LineWidth',1.5);
grid on; xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Reference','PID'); title('Case A: no acceleration limit');

figure; plot(t,caseA.a,'b','LineWidth',1.5);
grid on; xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
title('Case A: commanded acceleration');

figure; plot(t,vref,'k--',t,caseB.v,'b','LineWidth',1.5);
grid on; xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Reference','PID'); title('Case B: acceleration limited');

figure; plot(t,caseB.a,'b','LineWidth',1.5);
grid on; xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
title('Case B: saturated acceleration');

function out = simulateSpeedPID(vref, Ts, gains, L, accelLim)
N = numel(vref);
x = [0; 0; 0; 0];           % [x; y; psi; v]
Iev = 0;
ePrev = 0;

out.x  = zeros(4,N);
out.v  = zeros(1,N);
out.a  = zeros(1,N);
out.ev = zeros(1,N);

aMin = accelLim(1);
aMax = accelLim(2);

for k = 1:N
    % Reference-rate feedforward.
    if k == 1
        dvref = 0;
    else
        dvref = (vref(k) - vref(k-1))/Ts;
    end

    % Error and derivative terms.
    e = vref(k) - x(4);
    if k == 1
        dev = 0;
    else
        dev = (e - ePrev)/Ts;
    end

    % PID law: unsaturated command with the current integrator I_k.
    aUnsat = gains.Kp*e + gains.Ki*Iev + ...
             gains.Kd*dev + gains.Kff*dvref;
    a = min(max(aUnsat, aMin), aMax);

    % Conditional anti-windup: I_{k+1} = I_k + Ts*e unless the command
    % is saturated and the error would drive it further into saturation.
    holdInt = (aUnsat >= aMax && e >= 0) || ...
              (aUnsat <= aMin && e <= 0);
    if ~holdInt
        Iev = Iev + Ts*e;
    end

    % Straight-line kinematic update with delta_k = 0.
    x = x + Ts*[x(4)*cos(x(3));
                x(4)*sin(x(3));
                (x(4)/L)*tan(0);
                a];

    out.x(:,k) = x;
    out.v(k) = x(4);
    out.a(k) = a;
    out.ev(k) = e;
    ePrev = e;
end
end
