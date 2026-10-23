% main.m - 4-DOF Robot Kinematic Analysis
clear; clc; close all;

% 1. ROBOT PARAMETERS
L1 = 0.5; L3 = 0.4; L4 = 0.2;
L_dim = [L1, L3, L4];
alpha = pi/6; % 30 degrees inclination

% 2. TRAJECTORY GENERATION (Sinusoidal Surface)
t = 0:0.1:10; % Time vector
% (Insert surface parameterization here later: S, Sp, Spp arrays)

% Dummy target for initialization test
S_target = [0.5; 0.3; 0.6; 0];
Sp_target = [0.1; 0; 0; 0];
Spp_target = [0; 0; 0; 0];

% 3. SOLVE INVERSE KINEMATICS
[Q_calc, Qp_calc, Qpp_calc] = InverseKinematics(S_target, Sp_target, Spp_target, L_dim, alpha);

% 4. VERIFY WITH FORWARD KINEMATICS
[M04, W04, H04] = ForwardKinematics(Q_calc, Qp_calc, Qpp_calc, L_dim, alpha);

% Display results
disp('Calculated Tool Position (from M04):');
disp(M04(1:3, 4));

% 5. PLOT
PlotRobot4DOF(Q_calc, L_dim, alpha, 1, 'b', '-', 'Robot Kinematic State');