function [M04, W04, H04] = ForwardKinematics(Q, Qp, Qpp, L_dim, alpha)
% ForwardKinematics: Calculates Tool Pose, Velocity, and Acceleration matrices

% Retrieve M and L matrices from the base function
[M, L_global] = ComputeMatrices(Q, L_dim, alpha);
M04 = M{4}; % Final TCP pose

% Extract global operator matrices
L1_0 = L_global{1}; L2_0 = L_global{2}; 
L3_0 = L_global{3}; L4_0 = L_global{4};

q1p = Qp(1); q2p = Qp(2); q3p = Qp(3); q4p = Qp(4);
q1pp = Qpp(1); q2pp = Qpp(2); q3pp = Qpp(3); q4pp = Qpp(4);

% --- VELOCITY MATRIX (W) ---
W04 = L1_0*q1p + L2_0*q2p + L3_0*q3p + L4_0*q4p;

% --- ACCELERATION MATRIX (H) ---
H_norm = (L1_0^2)*(q1p^2) + (L2_0^2)*(q2p^2) + (L3_0^2)*(q3p^2) + (L4_0^2)*(q4p^2);
H_tang = L1_0*q1pp + L2_0*q2pp + L3_0*q3pp + L4_0*q4pp;

H_coriolis = 2 * (L1_0 * L2_0 * q1p * q2p + ...
    L1_0 * L3_0 * q1p * q3p + ...
    L1_0 * L4_0 * q1p * q4p + ...
    L2_0 * L3_0 * q2p * q3p + ...
    L2_0 * L4_0 * q2p * q4p + ...
    L3_0 * L4_0 * q3p * q4p);

H04 = H_norm + H_tang + H_coriolis;
end