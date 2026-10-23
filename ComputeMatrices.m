function [M, L_global] = ComputeMatrices(Q, L_dim, alpha)
% ComputeMatrices: Generates Pose (M) and Operator (L) matrices

q1 = Q(1); q2 = Q(2); q3 = Q(3); q4 = Q(4);
L1 = L_dim(1); L3 = L_dim(2); L4 = L_dim(3);

% --- 1. LOCAL OPERATOR MATRICES (L) ---
L1_loc = [0 -1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0]; % Revolute Z
L2_loc = [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0];  % Prismatic X
L3_loc = [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0]; % Revolute Y
L4_loc = [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0]; % Revolute Y

% --- 2. RELATIVE POSE MATRICES (M) ---
M01 = [cos(q1) -sin(q1) 0 0; sin(q1) cos(q1) 0 0; 0 0 1 L1; 0 0 0 1];
M12 = [cos(alpha) 0 -sin(alpha) q2*cos(alpha); 0 1 0 0; sin(alpha) 0 cos(alpha) q2*sin(alpha); 0 0 0 1];
M23 = [cos(q3) 0 -sin(q3) L3*cos(q3); 0 1 0 0; sin(q3) 0 cos(q3) L3*sin(q3); 0 0 0 1];
M34 = [cos(q4) 0 -sin(q4) L4*cos(q4); 0 1 0 0; sin(q4) 0 cos(q4) L4*sin(q4); 0 0 0 1];

% --- 3. GLOBAL POSE MATRICES ---
M02 = M01 * M12;
M03 = M02 * M23;
M04 = M03 * M34;

% --- 4. GLOBAL OPERATOR MATRICES ---
% Transforming L into the base frame: L_i(0) = M_0i * L_i(loc) * inv(M_0i)
L1_0 = L1_loc; 
L2_0 = M01 * L2_loc * inv(M01);
L3_0 = M02 * L3_loc * inv(M02);
L4_0 = M03 * L4_loc * inv(M03);

% Pack outputs into cell arrays for easy access
M = {M01, M02, M03, M04};
L_global = {L1_0, L2_0, L3_0, L4_0};
end