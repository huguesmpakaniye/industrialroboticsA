function Fq = InverseDynamics(Q, Qp, Qpp, L_dim, alpha, Masses)
% InverseDynamics: motor efforts Fq = [tau1; F2; tau3; tau4] with the extended Jacobian
% method of Prof. Giberti's course:
%     Fq = (Je'*M*Je)*Qpp + (Je'*M*Je_dot)*Qp - Je'*Fse
% Point masses placed at the REAL centres of mass of the CAD parts (STEP files),
% i.e. at the same places as in the Simscape model:
%   body 1 (q1): Part2 on the q1 axis + the 2 plates Part3
%                (123.5 mm along the prismatic axis from the pivot, +-135 mm aside)
%   body 2 (q2): Part4 on the prismatic axis, 126.2 mm BEFORE the elbow (moves with the whole q2)
%                + motor/reducer of q3 on the elbow axis            (Masses.mot3, optional)
%   body 3 (q3): Part5 on L3, 88.8 mm from the elbow
%                + motor/reducer of q4 at mid L3                    (Masses.mot4, optional)
%   body 4 (q4): Part6 on the tool axis, 47.2 mm from the wrist
% Masses.m1..m4 = masses of the 4 bodies (without motors), [kg]

% --- 1. CENTRES OF GRAVITY FROM THE CAD -----------------------------------------
d1u = 0.1235;                          % plates CoG along the prismatic axis from the pivot [m]
d1l = 0.135;                           % plates CoG lateral offset (each side) [m]
k1  = 2*2.174/(12.915 + 2*2.174);      % mass fraction of the 2 plates in body 1 (CAD volumes)
d2  = -0.1262;                         % Part4 CoG from the elbow, along the prismatic axis [m]
d3  = 0.0888;                          % Part5 CoG from the elbow, along L3 [m]
dm4 = L_dim(2)/2;                      % motor of q4 at mid L3 [m] (as in Simscape)
d4  = 0.0472;                          % Part6 CoG from the wrist, along L4 [m]

mot3 = 0;  mot4 = 0;
if isfield(Masses, 'mot3'), mot3 = Masses.mot3; end
if isfield(Masses, 'mot4'), mot4 = Masses.mot4; end
mPl = k1*Masses.m1/2;                  % each plate
mP2 = (1 - k1)*Masses.m1;              % Part2 (on the q1 axis)
m   = [mP2, mPl, mPl, Masses.m2, mot3, Masses.m3, mot4, Masses.m4];   % 8 point masses

% extended position vector Se(Q): [x; y; z] of every point mass (24 x 1)
calc_Se = @(q) pointsPos(q, L_dim(1), L_dim(2), alpha, d1u, d1l, d2, d3, dm4, d4);

% --- 2. EXTENDED JACOBIAN Je = dSe/dQ AND Je_dot (numerical, central differences) ---
Je     = jacNum(calc_Se, Q, 1e-6);
h      = 1e-5;                         % Je_dot*Qp = derivative of Je along Qp
Je_dot = (jacNum(calc_Se, Q + Qp*h, 1e-6) - jacNum(calc_Se, Q - Qp*h, 1e-6))/(2*h);

% --- 3. MASS MATRIX AND GRAVITY --------------------------------------------------
g   = 9.81;
M   = diag(kron(m, [1 1 1]));          % [m1 m1 m1 m2 m2 m2 ...]
Fse = zeros(3*numel(m), 1);
Fse(3:3:end) = -m*g;                   % gravity on the z components

% --- 4. INVERSE DYNAMICS -----------------------------------------------------------
Fq = (Je'*M*Je)*Qpp + (Je'*M*Je_dot)*Qp - Je'*Fse;
end

function Se = pointsPos(q, L1, L3, a, d1u, d1l, d2, d3, dm4, d4)
% [x; y; z] of the 8 point masses. In the arm plane: (radial, z); lateral = across the plane
u  = [cos(a); sin(a)];                         % prismatic axis
e3 = [cos(a + q(3)); sin(a + q(3))];           % L3 direction
e4 = [cos(a + q(3) + q(4)); sin(a + q(3) + q(4))];   % tool direction
P0 = [0; L1];                                  % pivot = axes intersection
E  = P0 + q(2)*u;                              % elbow
W  = E + L3*e3;                                % wrist
Pp = P0 + d1u*u;                               % plates CoG (in plane)
pts = [ 0,     0,    0;                        % Part2 (on the q1 axis, height irrelevant)
    Pp(1), +d1l, Pp(2);                    % plate 1
    Pp(1), -d1l, Pp(2);                    % plate 2
    [E + d2*u;   0].' * [1 0 0; 0 0 1; 0 1 0];    % Part4
    [E;          0].' * [1 0 0; 0 0 1; 0 1 0];    % motor q3
    [E + d3*e3;  0].' * [1 0 0; 0 0 1; 0 1 0];    % Part5
    [E + dm4*e3; 0].' * [1 0 0; 0 0 1; 0 1 0];    % motor q4
    [W + d4*e4;  0].' * [1 0 0; 0 0 1; 0 1 0] ];  % Part6
c1 = cos(q(1));  s1 = sin(q(1));
x  = pts(:,1)*c1 - pts(:,2)*s1;                % rotation q1 about z
y  = pts(:,1)*s1 + pts(:,2)*c1;
Se = reshape([x, y, pts(:,3)].', [], 1);
end

function J = jacNum(f, q, d)
% central-difference Jacobian of f at q
f0 = f(q);  J = zeros(numel(f0), numel(q));
for i = 1:numel(q)
    e = zeros(size(q));  e(i) = d;
    J(:,i) = (f(q + e) - f(q - e))/(2*d);
end
end