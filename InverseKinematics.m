function [Q, Qp, Qpp] = InverseKinematics(S, Sp, Spp, L_dim, alpha)
% InverseKinematics: Calculates Joint Position, Velocity, and Acceleration
% S, Sp, Spp are 4x1 vectors representing [X; Y; Z; Phi] and their derivatives

% 1. Position (using your existing geometric analytical solver)
Q = Robot4DOF_inv(S, L_dim, alpha); 

% 2. Extract Jacobian J and its time derivative J_dot
J = Robot4DOF_jac(Q, L_dim, alpha);

% 3. Velocity
Qp = J \ Sp; % Solves Qp = inv(J) * Sp

% 4. Acceleration
Jp = Robot4DOF_jacDot(Q, Qp, L_dim, alpha);
Qpp = J \ (Spp - Jp * Qp); % Solves Qpp = inv(J) * (Spp - Jp*Qp)
end