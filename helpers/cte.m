function e = cte(next,pass,p)
%next=next waypoint, pass = passed waypoint, p=current car point

%Linear line between the 2 waypoints
c(1,1:2) = polyfit([next(1,1) pass(1,1)],[next(2,1) pass(2,1)],1); %c(1)x-y+c(2)=0

%Check if p is on the left or right
if norm(p-next) < 1e-9   % vehicle is exactly at next waypoint
    e = 0;
    return;
end
a = (p-next)/norm(p-next); %Unit vector from next to p
b = (pass-next)/norm(pass-next); %Unit vector from next to pass

%Find angle of the two vectors (clamp dot product to [-1,1] to avoid acos(NaN))
theta = rad2deg(acos(min(1, max(-1, a'*b))));

%value=NEGATIVE--> point on the right of b -->CTE=POSITIVE
value = a(1)*b(2)-a(2)*b(1);

if abs(theta) == 0  %a and b are on the same direction
    e = 0;
    %disp('point on line');
elseif abs(theta) > 0 & value >= 0
    %disp('point on left');
    e = -abs(c(1,1)*p(1,1)-p(2,1)+c(1,2))/sqrt(c(1,1)^2+(-1)^2); %negative error
elseif abs(theta) > 0 & value < 0
    %disp('point on right');
    e = +abs(c(1,1)*p(1,1)-p(2,1)+c(1,2))/sqrt(c(1,1)^2+(-1)^2); %positive error
end