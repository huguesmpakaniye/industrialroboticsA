function Jp = Robot4DOF_jacDot(Q, Qp, L_dim, alpha)
% Robot4DOF_jacDot: Time derivative of the Jacobian matrix

q1 = Q(1); q2 = Q(2); q3 = Q(3); q4 = Q(4);
q1p = Qp(1); q2p = Qp(2); q3p = Qp(3); q4p = Qp(4);

L3 = L_dim(2); L4 = L_dim(3);

% Angles and radial extensions
phi = alpha + q3 + q4;
R = q2 * cos(alpha) + L3 * cos(alpha + q3) + L4 * cos(phi);

% Time derivatives of angles and extensions
phip = q3p + q4p;
Rp = q2p * cos(alpha) - L3 * sin(alpha + q3) * q3p - L4 * sin(phi) * phip;

% Initialize 4x4 matrix
Jp = zeros(4,4);

% Derivative of Row 1 (Xp terms)
Jp(1,1) = -Rp * sin(q1) - R * cos(q1) * q1p;
Jp(1,2) = -cos(alpha) * sin(q1) * q1p;
term1_3 = -L3 * cos(alpha + q3) * q3p - L4 * cos(phi) * phip;
base1_3 = -L3 * sin(alpha + q3) - L4 * sin(phi);
Jp(1,3) = term1_3 * cos(q1) - base1_3 * sin(q1) * q1p;
Jp(1,4) = -L4 * cos(phi) * phip * cos(q1) - L4 * sin(phi) * sin(q1) * q1p;

% Derivative of Row 2 (Yp terms)
Jp(2,1) = Rp * cos(q1) - R * sin(q1) * q1p;
Jp(2,2) = cos(alpha) * cos(q1) * q1p;
Jp(2,3) = term1_3 * sin(q1) + base1_3 * cos(q1) * q1p;
Jp(2,4) = -L4 * cos(phi) * phip * sin(q1) + L4 * sin(phi) * cos(q1) * q1p;

% Derivative of Row 3 (Zp terms)
Jp(3,1) = 0;
Jp(3,2) = 0;
Jp(3,3) = -L3 * sin(alpha + q3) * q3p - L4 * sin(phi) * phip;
Jp(3,4) = -L4 * sin(phi) * phip;

% Derivative of Row 4 (Phi terms)
Jp(4,1) = 0;
Jp(4,2) = 0;
Jp(4,3) = 0;
Jp(4,4) = 0;
end