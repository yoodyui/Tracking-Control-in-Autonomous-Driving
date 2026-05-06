function [W, velocities, headings] = waypoints_with_velocity_and_heading(n_points)
    % Generate parametric waypoints for a smooth soccer field-like shape

    % Segment 1: Straight line along the bottom side of the field
    X_bottom = linspace(30, 115, n_points/4);  % Straight line from left to right
    Y_bottom = linspace(0, 0, n_points/4);     % Y constant at 0
    
    % Segment 2: Right-side curve (smooth arc)
    r = 55;  % Radius for the curves
    theta_right = linspace(-pi/2, pi/2, n_points/4); 
    X_curve_right = r * cos(theta_right) + 115.02;  % X-coordinates of the right curve (shifted)
    Y_curve_right = r * sin(theta_right) + 55;      % Y-coordinates of the right curve (shifted)
    
    % Segment 3: Straight line along the top side of the field
    X_top = linspace(115.01, 30, n_points/4);  % Straight line from right to left
    Y_top = linspace(110, 110, n_points/4);    % Y constant at 110
    
    % Segment 4: Left-side curve (smooth arc)
    theta_left = linspace(pi/2, -pi/2, n_points/4);  
    X_curve_left = -r * cos(theta_left) + 29.9;  % X-coordinates of the left curve (shifted)
    Y_curve_left = r * sin(theta_left) + 55;     % Y-coordinates of the left curve (shifted)
    
    % Combine all segments into waypoints
    x = [X_bottom, X_curve_right, X_top, X_curve_left];
    y = [Y_bottom, Y_curve_right, Y_top, Y_curve_left];
    
    % Initialize velocities and headings
    velocities = zeros(1, n_points);
    headings = zeros(1, n_points); % Heading angles (radians)
    
    % Define speed limits
    max_speed = 50; % Maximum speed (m/s) (60m/s==216km/hr)
    min_speed = 1;  % Minimum speed (m/s)

    % Generate speed profile and headings
    for i = 1:n_points
        % Calculate heading using the past and next waypoints
        if i == 1
            % First waypoint: use last and second waypoints
            dx = x(i+1) - x(end);
            dy = y(i+1) - y(end);
        elseif i == n_points
            % Last waypoint: use second-to-last and first waypoints
            dx = x(1) - x(i-1);
            dy = y(1) - y(i-1);
        else
            % Middle waypoints: use previous and next waypoints
            dx = x(i+1) - x(i-1);
            dy = y(i+1) - y(i-1);
        end
        
        % Calculate heading angle (psi_ref)
        headings(i) = atan2(dy, dx);

        % Estimate curvature using three consecutive points
        if i == 1
            % Use the first and second points for the initial curvature
            p1 = [x(end), y(end)];
            p2 = [x(i), y(i)];
            p3 = [x(i+1), y(i+1)];
        elseif i == n_points
            % Use the last and second-to-last points for the final curvature
            p1 = [x(i-1), y(i-1)];
            p2 = [x(i), y(i)];
            p3 = [x(1), y(1)];
        else
            % Use consecutive points in the middle
            p1 = [x(i-1), y(i-1)];
            p2 = [x(i), y(i)];
            p3 = [x(i+1), y(i+1)];
        end

        % Calculate curvature
        curvature(i) = calculate_curvature(p1, p2, p3);

        % Assign speed dynamically based on curvature
        velocities(i) = max_speed - (max_speed - min_speed) * (min(curvature(i) / 0.03, 1));

    end
    

    % Combine waypoints, velocities, and headings
    W = [x; y];

end
