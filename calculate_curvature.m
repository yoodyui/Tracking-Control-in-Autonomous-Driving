function curvature = calculate_curvature(p1, p2, p3)
    % Calculate curvature from three points
    % p1, p2, p3 are [x, y] coordinates of three consecutive points

    % Triangle edge lengths
    a = norm(p2 - p3);
    b = norm(p1 - p3);
    c = norm(p1 - p2);

    % Semi-perimeter
    s = (a + b + c) / 2;

    % Area of the triangle using Heron's formula
    A = sqrt(s * (s - a) * (s - b) * (s - c));

    if A == 0
        curvature = 0; % Points are collinear
    else
        % Curvature is inverse of the circumradius
        curvature = 4 * A / (a * b * c);
    end
end