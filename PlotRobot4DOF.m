function PlotRobot4DOF(Q, L, alpha, fig_num, line_color, line_style, plot_title)
% PlotRobot4DOF: 3D Visualization of the 4 DOF robot

% Extract variables
q1 = Q(1); q2 = Q(2); q3 = Q(3); q4 = Q(4);
L1 = L(1); L3 = L(2); L4 = L(3);

% Define joints coordinates
x0 = 0; y0 = 0; z0 = 0;
x1 = 0; y1 = 0; z1 = L1;

R2 = q2 * cos(alpha);
z2 = L1 + q2 * sin(alpha);
x2 = R2 * cos(q1);
y2 = R2 * sin(q1);

R3 = R2 + L3 * cos(alpha + q3);
z3 = z2 + L3 * sin(alpha + q3);
x3 = R3 * cos(q1);
y3 = R3 * sin(q1);

phi = alpha + q3 + q4;
R4 = R3 + L4 * cos(phi);
z4 = z3 + L4 * sin(phi);
x4 = R4 * cos(q1);
y4 = R4 * sin(q1);

% Create or call the specified figure number[cite: 84]
figure(fig_num); 
hold on; grid on;

X = [x0, x1, x2, x3, x4];
Y = [y0, y1, y2, y3, y4];
Z = [z0, z1, z2, z3, z4];

% Plot the robot using the specified color and style
plot3(X, Y, Z, 'LineStyle', line_style, 'Marker', 'o', ...
    'LineWidth', 3, 'MarkerSize', 8, 'Color', line_color, ...
    'MarkerFaceColor', line_color);

% Draw the prismatic axis dashed line
R_axis = (q2 + 0.1) * cos(alpha);
Z_axis = L1 + (q2 + 0.1) * sin(alpha);
plot3([x1, R_axis*cos(q1)], [y1, R_axis*sin(q1)], [z1, Z_axis], ...
    '--k', 'LineWidth', 1);

axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(plot_title); 
view(3);

max_reach = L1 + q2 + L3 + L4;
xlim([-max_reach max_reach]);
ylim([-max_reach max_reach]);
zlim([0 max_reach]);
hold off;
end