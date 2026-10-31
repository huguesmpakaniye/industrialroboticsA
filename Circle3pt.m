function [P, R] = Circle3pt(P1, P2, P3, n)
% Circle3pt  n points (3 x n) on the circular arc P1 -> P2 -> P3 (3D, columns).
P1 = P1(:); P2 = P2(:); P3 = P3(:);
a = P2 - P1;  b = P3 - P1;
w = cross(a, b);
ctr = P1 + (cross(w, a)*(b.'*b) + cross(b, w)*(a.'*a)) / (2*(w.'*w));
R = norm(P1 - ctr);
u = (P1 - ctr)/R;
nv = w/norm(w);
v = cross(nv, u);
ang = @(Pq) mod(atan2((Pq-ctr).'*v, (Pq-ctr).'*u), 2*pi);
a2 = ang(P2);  a3 = ang(P3);
if a3 < a2                      % the arc must pass through P2
    th = linspace(0, a3 - 2*pi, n);
else
    th = linspace(0, a3, n);
end
P = zeros(3, n);
for i = 1:n
    P(:,i) = ctr + R*(cos(th(i))*u + sin(th(i))*v);
end
end

