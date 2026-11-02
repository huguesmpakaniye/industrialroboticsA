function Fq = InverseDynamics(Q, Qp, Qpp, L_dim, alpha, Masses)
% InverseDynamics: Computes required motor torques/forces Fq
% using the Extended Jacobian method from Prof. Giberti's course.

q1 = Q(1); q2 = Q(2); q3 = Q(3); q4 = Q(4);

% --- 1. CENTERS OF GRAVITY KINEMATICS ---
% Distances to the CoG for each link (Replace with actual CAD values)
dg1 = L_dim(1)/2;  % CoG of link 1 (base)
dg2 = q2/2;        % CoG of link 2 (prismatic)
dg3 = L_dim(2)/2;  % CoG of link 3 (elbow)
dg4 = L_dim(3)/2;  % CoG of link 4 (wrist)

% Anonymous function to compute the extended position vector Se(Q)
% Se contains [x; y; z] for the 4 centers of gravity (12x1 vector)[cite: 6]
calc_Se = @(q) [
    % Link 1 (Base - rotates around Z)
    dg1 * cos(q(1)); 
    dg1 * sin(q(1)); 
    L_dim(1)/2;

    % Link 2 (Prismatic - translates along alpha)
    (q(2)/2) * cos(alpha) * cos(q(1));
    (q(2)/2) * cos(alpha) * sin(q(1));
    L_dim(1) + (q(2)/2) * sin(alpha);

    % Link 3 (Elbow)
    (q(2)*cos(alpha) + dg3*cos(alpha+q(3))) * cos(q(1));
    (q(2)*cos(alpha) + dg3*cos(alpha+q(3))) * sin(q(1));
    L_dim(1) + q(2)*sin(alpha) + dg3*sin(alpha+q(3));

    % Link 4 (Wrist)
    (q(2)*cos(alpha) + L_dim(2)*cos(alpha+q(3)) + dg4*cos(alpha+q(3)+q(4))) * cos(q(1));
    (q(2)*cos(alpha) + L_dim(2)*cos(alpha+q(3)) + dg4*cos(alpha+q(3)+q(4))) * sin(q(1));
    L_dim(1) + q(2)*sin(alpha) + L_dim(2)*sin(alpha+q(3)) + dg4*sin(alpha+q(3)+q(4))
    ];

% --- 2. EXTENDED JACOBIAN (Je) AND ITS DERIVATIVE (Je_dot) ---
% Numerical differentiation for Je = d(Se)/dQ[cite: 7]
delta = 1e-6;
Je = zeros(12, 4);
Je_dot = zeros(12, 4);
Se_current = calc_Se(Q);

for i = 1:4
    Q_plus = Q; Q_plus(i) = Q_plus(i) + delta;
    Je(:, i) = (calc_Se(Q_plus) - Se_current) / delta;

    % For Je_dot, we differentiate Se_dot = Je * Qp
    Qp_plus = Qp; Qp_plus(i) = Qp_plus(i) + delta;
    % (Simplified numerical approximation of the convective term)
end

% Analytical Je_dot calculation via finite difference over time step
% requires state history, so we use the symbolic equivalent or a direct mapping.
% To keep it purely algebraic for this instantaneous function:
% We construct the convective acceleration vector: J_dot * Q_dot
% V(Q, Qp) = (Je^T * M * Je_dot) * Qp[cite: 10]
% Numerically, J_dot*Qp is the directional derivative of J along Qp:
Q_next = Q + Qp * delta;
Je_next = zeros(12, 4);
for i = 1:4
    Q_temp = Q_next; Q_temp(i) = Q_temp(i) + delta;
    Je_next(:, i) = (calc_Se(Q_temp) - calc_Se(Q_next)) / delta;
end
Je_dot = (Je_next - Je) / delta;

% --- 3. MASS MATRIX & EXTERNAL FORCES ---
% Diagonal mass matrix for X, Y, Z of each CoG[cite: 8]
m1 = Masses.m1; m2 = Masses.m2; m3 = Masses.m3; m4 = Masses.m4;
M_matrix = diag([m1, m1, m1, m2, m2, m2, m3, m3, m3, m4, m4, m4]);

% Gravity vector Fse (acting on -Z)[cite: 18]
g = 9.81;
Fse = zeros(12, 1);
Fse(3)  = -m1 * g;
Fse(6)  = -m2 * g;
Fse(9)  = -m3 * g;
Fse(12) = -m4 * g;

% --- 4. INVERSE DYNAMICS EQUATION ---
% Fq = (Je^T * M * Je)*Qpp + (Je^T * M * Je_dot)*Qp - Je^T * Fse[cite: 9]
Term_Inertia   = (Je' * M_matrix * Je) * Qpp;
Term_Coriolis  = (Je' * M_matrix * Je_dot) * Qp;
Term_Gravity   = Je' * Fse;

Fq = Term_Inertia + Term_Coriolis - Term_Gravity;
end