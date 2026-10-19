function S = Robot4DOF_dir(Q, L, alpha)
% Robot4DOF_dir: Forward kinematics for the 4 DOF robot
% Inputs:
%   Q     - Joint space vector [q1; q2; q3; q4] 
%           (q1, q3, q4 in radians; q2 in meters)
%   L     - Link lengths vector [L1; L3; L4] in meters
%   alpha - Constant inclination angle of the prismatic joint (radians)
% Output:
%   S     - Task space vector [Xp; Yp; Zp; phi]

% Extract variables
q1 = Q(1); q2 = Q(2); q3 = Q(3); q4 = Q(4);
L1 = L(1); L3 = L(2); L4 = L(3);

% Global orientation of the end-effector in the arm's plane (pitch)
phi = alpha + q3 + q4;

% Radial extension from the vertical base axis
R = q2 * cos(alpha) + L3 * cos(alpha + q3) + L4 * cos(phi);

% Initialize output vector
S = zeros(4,1);

% Calculate Cartesian coordinates[cite: 78]
S(1) = R * cos(q1);                                        % Xp
S(2) = R * sin(q1);                                        % Yp
S(3) = L1 + q2 * sin(alpha) + L3 * sin(alpha + q3) + L4 * sin(phi); % Zp
S(4) = phi;                                                % Orientation phi
end