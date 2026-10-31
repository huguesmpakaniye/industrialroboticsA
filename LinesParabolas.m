function [pos, vel, acc, ok] = LinesParabolas(q, tau, accMax, t)
% LinesParabolas  Linear segments with parabolic blends for ONE joint.
%   q      : n nodes (n >= 4), joint values at the nodes
%   tau    : n-1 segment durations (same for all joints!)
%   accMax : blend acceleration magnitude of this joint (> 0)
%   t      : evaluation times (0 ... sum(tau))
% The first and last node are reached exactly (start/end at rest).
% Interior nodes are NOT reached: the parabola cuts the corner.
% ok = false if a blend does not fit (raise accMax, slow down, or add nodes).

q = q(:); tau = tau(:); n = numel(q);
pos = []; vel = []; acc = []; ok = true;

v  = zeros(n-1,1);      % velocity of each linear segment
tb = zeros(n,1);        % blend duration at each node
a  = zeros(n,1);       % blend acceleration at each node

% --- first blend (starts at rest at node 1, ends inside segment 1)
dq = q(2) - q(1);
if dq ~= 0
    a(1) = sign(dq)*accMax;
    rad  = tau(1)^2 - 2*dq/a(1);
    if rad < 0, ok = false; return; end
    tb(1) = tau(1) - sqrt(rad);
end
v(1) = dq/(tau(1) - tb(1)/2);

% --- last blend (ends at rest at node n)
dq = q(n) - q(n-1);
if dq ~= 0
    a(n) = -sign(dq)*accMax;
    rad  = tau(n-1)^2 + 2*dq/a(n);
    if rad < 0, ok = false; return; end
    tb(n) = tau(n-1) - sqrt(rad);
end
v(n-1) = dq/(tau(n-1) - tb(n)/2);

% --- interior segment velocities and blends (centred on the node time)
for j = 2:n-2
    v(j) = (q(j+1) - q(j))/tau(j);
end
for j = 2:n-1
    dv = v(j) - v(j-1);
    if dv ~= 0
        a(j)  = sign(dv)*accMax;
        tb(j) = dv/a(j);
    end
end

% --- duration of the linear parts
tl = zeros(n-1,1);
tl(1)   = tau(1)   - tb(1)   - tb(2)/2;
tl(n-1) = tau(n-1) - tb(n-1)/2 - tb(n);
for j = 2:n-2
    tl(j) = tau(j) - tb(j)/2 - tb(j+1)/2;
end
if min(tl) < -1e-12, ok = false; return; end
tl = max(tl, 0);

% --- chain of intervals: blend1, line1, blend2, line2, ... blend n
m   = 2*n - 1;
dur = zeros(m,1);  aint = zeros(m,1);
dur(1:2:end)  = tb;  aint(1:2:end) = a;
dur(2:2:end)  = tl;
Tk = [0; cumsum(dur)];
P0 = zeros(m,1); V0 = zeros(m,1);
P0(1) = q(1); V0(1) = 0;
for k = 1:m-1
    P0(k+1) = P0(k) + V0(k)*dur(k) + 0.5*aint(k)*dur(k)^2;
    V0(k+1) = V0(k) + aint(k)*dur(k);
end

% --- evaluate
t = t(:);
k = sum(bsxfun(@ge, t, Tk(1:end-1).'), 2);   % last interval whose start <= t
k = min(max(k,1), m);
dt  = t - Tk(k);
pos = P0(k) + V0(k).*dt + 0.5*aint(k).*dt.^2;
vel = V0(k) + aint(k).*dt;
acc = aint(k);
end