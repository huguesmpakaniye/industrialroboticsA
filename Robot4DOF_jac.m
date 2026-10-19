function J = Robot4DOF_jac(Q, L, alpha)
% Robot4DOF_jac: Jacobian matrix for the 4 DOF robot[cite: 78]
% Inputs:
%   Q     - Joint space vector [q1; q2; q3; q4]
%   L     - Link lengths vector [L1; L3; L4]
%   alpha - Constant inclination angle of the prismatic joint
% Output:
%   J     - 4x4 Jacobian matrix

q1 = Q(1); q2 = Q(2); q3 = Q(3); q4 = Q(4);
L3 = L(2); L4 = L(3);

phi = alpha + q3 + q4;
R = q2 * cos(alpha) + L3 * cos(alpha + q3) + L4 * cos(phi);

% Initialize 4x4 Jacobian matrix
J = zeros(4,4);

% Row 1: Derivatives of Xp w.r.t q1, q2, q3, q4[cite: 78]
J(1,1) = -R * sin(q1);
J(1,2) = cos(alpha) * cos(q1);
J(1,3) = (-L3 * sin(alpha + q3) - L4 * sin(phi)) * cos(q1);
J(1,4) = -L4 * sin(phi) * cos(q1);

% Row 2: Derivatives of Yp w.r.t q1, q2, q3, q4
J(2,1) = R * cos(q1);
J(2,2) = cos(alpha) * sin(q1);
J(2,3) = (-L3 * sin(alpha + q3) - L4 * sin(phi)) * sin(q1);
J(2,4) = -L4 * sin(phi) * sin(q1);

% Row 3: Derivatives of Zp w.r.t q1, q2, q3, q4
J(3,1) = 0;
J(3,2) = sin(alpha);
J(3,3) = L3 * cos(alpha + q3) + L4 * cos(phi);
J(3,4) = L4 * cos(phi);

% Row 4: Derivatives of global orientation phi w.r.t q1, q2, q3, q4
J(4,1) = 0;
J(4,2) = 0;
J(4,3) = 1;
J(4,4) = 1;
end