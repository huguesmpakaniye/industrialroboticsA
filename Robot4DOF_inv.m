function Q = Robot4DOF_inv(S, L, alpha)
% Robot4DOF_inv: Inverse kinematics for the 4 DOF robot[cite: 80]
% Inputs:
%   S     - Task space target vector [Xp; Yp; Zp; phi]
%   L     - Link lengths vector [L1; L3; L4]
%   alpha - Constant inclination angle of the prismatic joint
% Output:
%   Q     - Joint space vector [q1; q2; q3; q4]

% Extract variables
Xp = S(1); Yp = S(2); Zp = S(3); phi = S(4);
L1 = L(1); L3 = L(2); L4 = L(3);

% Initialize output vector
Q = zeros(4,1);

% Step 1: Base rotation q1
Q(1) = atan2(Yp, Xp); % atan2 handles all quadrants[cite: 76, 80]

% Step 2: Wrist center decoupling (Point W coordinates)[cite: 40, 41]
% R_w is the radial distance of the wrist center from the z-axis
R_w = sqrt(Xp^2 + Yp^2) - L4 * cos(phi);
Z_w = Zp - L4 * sin(phi);

% Step 3: Solve for elbow rotation q3
% Projection of the wrist coordinates to isolate q3
numerator = (Z_w - L1) * cos(alpha) - R_w * sin(alpha);
Q(3) = asin(numerator / L3); 

% Step 4: Solve for prismatic joint q2
% Projection along the inclined axis alpha
Q(2) = R_w * cos(alpha) + (Z_w - L1) * sin(alpha) - L3 * cos(Q(3));

% Step 5: Solve for wrist rotation q4
Q(4) = phi - alpha - Q(3);
end