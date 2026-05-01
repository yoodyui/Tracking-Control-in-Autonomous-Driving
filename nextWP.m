function [next,pass,next_index,pass_index] = nextWP(W,p)

for i = 1:length(W)
    dist(i) = dis(p,W(:,i));
end

%Find the nearest waypoint to current point
index_min = find(dist==min(dist));
index_min = index_min(1);   % take first if multiple equidistant waypoints

%Define the successive waypoint
index_next = index_min+1;
if index_next==length(W)+1
    index_next = 1;
end

%Define the preceding waypoint
index_past = index_min-1;
if index_past==0
    index_past = length(W);
end

%Check it W==p or not
if min(dist)==0 %W==p

    pass_index = index_min;
    next_index = index_next;

else %W~=p

    %Test to see two angles angle(PWcWn) and angle(WcWnP)
    %v_1 = vector from closet waypoint to current point
    v_1(1:2,1) = (p-W(:,index_min))/norm(p-W(:,index_min));
    %v_2 = vector from closet waypoint to next index waypoint
    v_2(1:2,1) = (W(:,index_next)-W(:,index_min))/norm(W(:,index_next)-W(:,index_min));
    %Find cross angle bt 2 vectors
    ang1 = rad2deg(acos(v_1(:,1)'*v_2(:,1)));

    th = 90; %Thredshold

    if ang1<th %p lies in between next and closest waypoints
        pass_index = index_min;
        next_index = index_next;
        % disp('case1')

    else  %p lies in between past and closest waypoints
        pass_index = index_past;
        next_index = index_min;
        % disp('case2')
    end

end

pass = W(:,pass_index);
next = W(:,next_index);


% % %
% figure; plot(W(1,index_min),W(2,index_min),'*'); hold on; plot(W(1,index_next),W(2,index_next),'ob');
% plot(W(1,index_past),W(2,index_past),'or');
% plot(p(1,1),p(2,1),'x');

%plot([W(1,next_index) W(1,pass_index) p(1,1) W(1,next_index)],[W(2,next_index) W(2,pass_index) p(2,1) W(2,next_index)]);
%hold on; quiver(W(1,index_min),W(2,index_min),v_2(1,1),v_2(2,1));
%hold on; quiver(W(1,index_min),W(2,index_min),v_1(1,1),v_1(2,1));