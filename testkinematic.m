% TEST SCRIPT: 4 DOF Robot Direct and Inverse Kinematics[cite: 81]
clc; clear all; close all;

% 1. Define geometric parameters 
L1 = 0.200; % Base height in meters
L3 = 0.250; % Forearm length in meters
L4 = 0.100; % Wrist/TCP length in meters
alpha = pi/6; % 30 degrees inclination (prismatic joint)
L = [L1; L3; L4];

% 2. Define a test joint configuration
% Q = [base_rot; prism_ext; elbow_rot; wrist_rot]
Q_test = [pi/4; 0.150; pi/6; -pi/12]; 

disp('--- ORIGINAL JOINT VECTOR (Q_test) ---');
disp(Q_test);

% 3. Test Forward Kinematics (MCD)[cite: 81]
S_calc = Robot4DOF_dir(Q_test, L, alpha);
disp('--- CALCULATED END-EFFECTOR POSE (S_calc) ---');
disp('Xp, Yp, Zp (meters), Phi (rad):');
disp(S_calc);

% 4. Test Inverse Kinematics (MCI)[cite: 81]
Q_calc = Robot4DOF_inv(S_calc, L, alpha);
disp('--- CALCULATED JOINT VECTOR (Q_calc) ---');
disp(Q_calc);

% 5. Verification 
disp('--- ERROR (Should be close to zero) ---');
error_vector = Q_test - Q_calc;
disp(error_vector);

% 6. PLOT BOTH ROBOTS IN SEPARATE FIGURES
% Plot the original robot (Forward Kinematics) on Figure 1
PlotRobot4DOF(Q_test, L, alpha, 1, 'b', '-', 'Robot Forward Kinematic (Original Configuration)');

% Plot the calculated robot (Inverse Kinematics) on Figure 2
PlotRobot4DOF(Q_calc, L, alpha, 2, 'r', '-', 'Robot Inverse Kinematic (Calculated Configuration)');