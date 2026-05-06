function [next, pass, next_index, pass_index] = nextForwardWP(W, p, i_prev, window)
% nextWP  Return the "next" and "pass" waypoints near position p.
% 
% Usage:
%   [next,pass,next_index,pass_index] = nextWP(W,p)              % original: nearest-based
%   [next,pass,next_index,pass_index] = nextWP(W,p,i_prev)       % forward-only (default window=25)
%   [next,pass,next_index,pass_index] = nextWP(W,p,i_prev,window)% forward-only with custom window
%
% Inputs:
%   W       2xN waypoint matrix [x;y]
%   p       2x1 current position
%   i_prev  (optional) last NEXT index returned in previous call
%   window  (optional) number of points to search ahead of i_prev
%
% Outputs:
%   next, pass       2x1 waypoints
%   next_index       index of NEXT waypoint (monotone if i_prev is given)
%   pass_index       index of PASS (previous) waypoint

N = size(W,2);

% ---------- Forward-only mode (preferred) ----------
if nargin >= 3 && ~isempty(i_prev)
    if nargin < 4 || isempty(window), window = 25; end

    % forward-only candidate set with wrap-around
    cand = mod((i_prev+1):(i_prev+window)-1, N) + 1;

    % choose the nearest among forward candidates
    d2 = sum( (W(:,cand) - p).^2, 1 );
    [~, j] = min(d2);
    next_index = cand(j);

    % previous index (one behind next), wrapped
    pass_index = mod(next_index-2, N) + 1;

    next = W(:, next_index);
    pass = W(:, pass_index);
    return;
end

% ---------- Original nearest-based behavior (backward compatible) ----------
% vectorized distance
dist = sqrt(sum((W - p).^2, 1));
[dist_min, index_min] = min(dist);

% indices around the nearest point
index_next = index_min + 1;  if index_next > N, index_next = 1; end
index_past = index_min - 1;  if index_past < 1, index_past = N; end

if dist_min == 0
    pass_index = index_min;
    next_index = index_next;
else
    % decide which side of the nearest point p is on
    v1 = p - W(:,index_min);         n1 = norm(v1); if n1>0, v1 = v1/n1; end
    v2 = W(:,index_next) - W(:,index_min); n2 = norm(v2); if n2>0, v2 = v2/n2; end
    ang1 = acosd( max(-1,min(1, v1.'*v2)) );   % robust dot→angle

    if ang1 < 90      % p lies between nearest and next
        pass_index = index_min;
        next_index = index_next;
    else               % p lies between past and nearest
        pass_index = index_past;
        next_index = index_min;
    end
end

next = W(:, next_index);
pass = W(:, pass_index);
end
