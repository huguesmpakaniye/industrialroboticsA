% projectA_main.m - Project A: letter "A" on the sinusoidal cylinder
% Stages: (1) circular approach  (2) trace the A  (3) linear return home
% Uses: Robot4DOF_inv, InverseKinematics, ForwardKinematics, Robot4DOF_jac,
%       cycloidal, PlotRobot4DOF   (+ new: LinesParabolas, Circle3pt)
% Units: m, rad, s.   All numbers marked [ASSUMPTION] are yours to change.
clear; clc; close all;

%% 0. SURFACE AND ROBOT --------------------------------------------------
R0 = 0.250;  Amp = 0.01125;  Hs = 0.200;  kr = 2.5;      % given (A, H, R0); kr ripples over H [ASSUMPTION]
L1 = 0.10;  L3 = 0.16;  L4 = 0.06;                        % [ASSUMPTION] L4 includes the pen
L_dim = [L1 L3 L4];  alpha = pi/6;                        % [ASSUMPTION] 30 deg prismatic axis
zb = L1;  zc = zb + Hs/2;  th0 = 0;                       % surface bottom, letter centre height, azimuth
rS   = @(z) R0 + Amp*sin(2*pi*kr*(z - zb)/Hs);            % surface radius r(z)
drS  = @(z) Amp*2*pi*kr/Hs*cos(2*pi*kr*(z - zb)/Hs);      % r'(z)
phiS = @(z) -atan(drS(z));                                % tool pitch = surface normal

S_home = [0.20; 0; zb + Hs/2; 0];                         % [ASSUMPTION] home: inside the wall, tool horizontal
Q_home = Robot4DOF_inv(S_home, L_dim, alpha);

%% 1. TIMING / LIMITS (placeholders until you pick motors + reducers) -----
vTool  = 0.010;                 % tool speed on the letter [m/s]  (the wrist q4 limits this, see sweep)
dsNode = 0.002;                 % node spacing along the letter [m]
accMax = [5; 1; 10; 20];        % blend accelerations [rad/s2; m/s2; rad/s2; rad/s2] -> replace by motor limits
Tcirc  = 2.0;  nCirc = 13;      % circular approach: duration, spline nodes
Tret   = 3.0;                   % return duration
dt     = 0.01;                  % output sample time

%% 2. LETTER NODES (A = 60 wide x 80 tall) ---------------------------------
% vertices in the letter plane [x;y] (m), origin = letter centre
%        LF      CL       CR      AP     RF
V = [-30  -18.75  18.75    0     30;
     -40  -10    -10      40    -40]*1e-3;
order = [1 2 3 4 2 3 5];        % LF-CL-CR-AP-CL-CR-RF: no 180 deg reversal, crossbar traced twice
x = []; y = [];
for i = 1:numel(order)-1
    a = V(:,order(i)); b = V(:,order(i+1));
    n = max(1, round(norm(b - a)/dsNode));
    s = (0:n-1)/n;
    x = [x, a(1) + (b(1)-a(1))*s];  %#ok<AGROW>
    y = [y, a(2) + (b(2)-a(2))*s];  %#ok<AGROW>
end
x(end+1) = V(1,order(end));  y(end+1) = V(2,order(end));
z  = zc + y;  th = th0 + x/R0;  r = rS(z);                % wrap onto the surface
S_let = [r.*cos(th); r.*sin(th); z; phiS(z)];             % 4 x N task-space nodes
N = size(S_let, 2);

Q_let = zeros(4, N);
for i = 1:N
    Q_let(:,i) = Robot4DOF_inv(S_let(:,i), L_dim, alpha);
end
assert(isreal(Q_let), 'Letter not reachable: change L3/L4/alpha/L1');

%% 3. STAGE 2 - TRACE THE A: lines and parabolas in joint space -----------
chord = sqrt(sum(diff(S_let(1:3,:), 1, 2).^2, 1));        % 1 x (N-1)
ramp  = ones(1, N-1);  ramp(1:3) = [3 2 1.4];  ramp(end-2:end) = [1.4 2 3];
okAll = false;
for vTry = vTool*[1 0.8 0.6 0.4]                          % slow down until every blend fits
    tau = (chord./vTry).*ramp;
    okAll = true;
    for j = 1:4
        [~,~,~,ok] = LinesParabolas(Q_let(j,:), tau, accMax(j), 0);
        okAll = okAll && ok;
    end
    if okAll, break; end
end
assert(okAll, 'No feasible blend: raise accMax or add auxiliary nodes');
fprintf('Letter: v_tool = %.1f mm/s, duration = %.1f s\n', vTry*1e3, sum(tau));
tL = (0:dt:sum(tau)).';
if tL(end) < sum(tau) - 1e-9, tL(end+1) = sum(tau); end   % hit the last node exactly
QL = zeros(4, numel(tL)); QpL = QL; QppL = QL;
for j = 1:4
    [QL(j,:), QpL(j,:), QppL(j,:)] = LinesParabolas(Q_let(j,:), tau, accMax(j), tL);
end

%% 4. STAGE 1 - CIRCULAR APPROACH: joint-space cubic splines --------------
P1 = S_home(1:3);  P3 = S_let(1:3,1);                     % home -> first point of the A
Mid = (P1 + P3)/2;  ch = P3 - P1;
up = [0;0;1];  perp = up - (up.'*ch)/(ch.'*ch)*ch;  perp = perp/norm(perp);
P2 = Mid + 0.010*perp;                                    % [ASSUMPTION] 10 mm bulge keeps the tip off the wall
[Parc, Rc] = Circle3pt(P1, P2, P3, nCirc);
S_arc = [Parc; linspace(S_home(4), S_let(4,1), nCirc)];
Q_arc = zeros(4, nCirc);
for i = 1:nCirc, Q_arc(:,i) = Robot4DOF_inv(S_arc(:,i), L_dim, alpha); end
tc  = linspace(0, Tcirc, nCirc);
pp  = spline(tc, [zeros(4,1), Q_arc, zeros(4,1)]);        % clamped: zero velocity at both ends
c   = pp.coefs;
ppd  = mkpp(pp.breaks, [3*c(:,1) 2*c(:,2) c(:,3)], 4);
ppdd = mkpp(pp.breaks, [6*c(:,1) 2*c(:,2)], 4);
tC  = (0:dt:Tcirc);
QC = ppval(pp, tC);  QpC = ppval(ppd, tC);  QppC = ppval(ppdd, tC);

%% 5. STAGE 3 - LINEAR RETURN: straight line in task space, cycloidal law --
S_end = S_let(:,end);  dS = S_home - S_end;
tR = 0:dt:Tret;  QR = zeros(4,numel(tR)); QpR = QR; QppR = QR;
for k = 1:numel(tR)
    m = cycloidal(tR(k), Tret, 0, 1);                     % path coordinate s: 0 -> 1
    [QR(:,k), QpR(:,k), QppR(:,k)] = InverseKinematics(S_end + m.pos*dS, m.vel*dS, m.acc*dS, L_dim, alpha);
end

%% 6. ASSEMBLE -------------------------------------------------------------
T   = [tC, Tcirc + tL(2:end).', Tcirc + tL(end) + tR(2:end)];
Q   = [QC,   QL(:,2:end),   QR(:,2:end)];
Qp  = [QpC,  QpL(:,2:end),  QpR(:,2:end)];
Qpp = [QppC, QppL(:,2:end), QppR(:,2:end)];

%% 7. VERIFY WITH FORWARD KINEMATICS --------------------------------------
M = numel(T);  tip = zeros(3,M);  vtip = zeros(3,M);  atip = zeros(3,M);
for k = 1:M
    [M04, W04, H04] = ForwardKinematics(Q(:,k), Qp(:,k), Qpp(:,k), L_dim, alpha);
    p = M04(:,4);  tip(:,k) = p(1:3);
    vv = W04*p;  vtip(:,k) = vv(1:3);
    aa = H04*p;  atip(:,k) = aa(1:3);
end
kmid = round(M/2);                                        % consistency: W04 vs J (needs the fixed ComputeMatrices)
J = Robot4DOF_jac(Q(:,kmid), L_dim, alpha);
[~, Wm] = ForwardKinematics(Q(:,kmid), Qp(:,kmid), Qpp(:,kmid), L_dim, alpha);
vm = Wm*[tip(:,kmid); 1];
fprintf('Velocity check |W*p - J*qp| = %.2e (must be ~0)\n', norm(vm(1:3) - J(1:3,:)*Qp(:,kmid)));

rho  = sqrt(tip(1,:).^2 + tip(2,:).^2);
gap  = rS(tip(3,:)) - rho;                                % >0: tip is inside the wall, 0: on the wall
iLet = (T >= Tcirc) & (T <= Tcirc + tL(end));
fprintf('On the letter: max |radial error| = %.2f mm\n', 1e3*max(abs(gap(iLet))));
fprintf('Min wall gap on arc/return (excluding last/first 0.2 s) = %.1f mm\n', ...
        1e3*min(gap(T>0.05 & T<Tcirc-0.2)));
fprintf('Peak joint speeds   : %s\n', num2str(max(abs(Qp),[],2).', '%8.3f'));
fprintf('Peak joint accel.   : %s\n', num2str(max(abs(Qpp),[],2).', '%8.3f'));
fprintf('Joint ranges        : q1 [%.2f %.2f]  q2 [%.3f %.3f]  q3 [%.2f %.2f]  q4 [%.2f %.2f]\n', ...
        min(Q(1,:)),max(Q(1,:)),min(Q(2,:)),max(Q(2,:)),min(Q(3,:)),max(Q(3,:)),min(Q(4,:)),max(Q(4,:)));

%% 8. PLOTS ----------------------------------------------------------------
figure(1); clf; hold on; grid on; axis equal;
[Zs, Ts] = meshgrid(linspace(zb, zb+Hs, 60), linspace(-0.6, 0.6, 60));
surf(rS(Zs).*cos(Ts), rS(Zs).*sin(Ts), Zs, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot3(tip(1,:), tip(2,:), tip(3,:), 'r', 'LineWidth', 1.5);
plot3(S_let(1,:), S_let(2,:), S_let(3,:), 'k.', 'MarkerSize', 4);
PlotRobot4DOF(Q_home, L_dim, alpha, 1, 'b', '-', 'Project A trajectory');
PlotRobot4DOF(Q_let(:,1), L_dim, alpha, 1, 'g', '-', 'Project A trajectory');
xlim([-0.05 0.3]); ylim([-0.2 0.2]); zlim([0 0.35]);
figure(2); clf;
sp = {'q_1 [rad]','q_2 [m]','q_3 [rad]','q_4 [rad]'};
for j = 1:4
    subplot(4,3,3*j-2); plot(T, Q(j,:));   grid on; ylabel(sp{j}); if j==1, title('position'); end
    subplot(4,3,3*j-1); plot(T, Qp(j,:));  grid on;               if j==1, title('velocity');  end
    subplot(4,3,3*j);   plot(T, Qpp(j,:)); grid on;               if j==1, title('acceleration'); end
end
xlabel('t [s]');

%% 9. DYNAMIC ANALYSIS (Inverse Dynamics) ----------------------------------
% Mass properties (from SolidWorks data)
Masses.m1 = 4.50;  % Base [kg]
Masses.m2 = 2.10;  % Prismatic link [kg]
Masses.m3 = 1.80;  % Elbow [kg]
Masses.m4 = 0.75;  % Wrist + tool [kg]

% Initialize motor efforts array
Fq_history = zeros(4, M);

for k = 1:M
    Fq_history(:, k) = InverseDynamics(Q(:,k), Qp(:,k), Qpp(:,k), L_dim, alpha, Masses);
end

% --- Check maximum efforts (Motor Sizing) ---
max_efforts = max(abs(Fq_history), [], 2);
fprintf('\n--- DYNAMIC RESULTS (Motor Sizing) ---\n');
fprintf('Max Tau 1 (Base)        : %.2f Nm\n', max_efforts(1));
fprintf('Max Force 2 (Prismatic) : %.2f N\n', max_efforts(2));
fprintf('Max Tau 3 (Elbow)       : %.2f Nm\n', max_efforts(3));
fprintf('Max Tau 4 (Wrist)       : %.2f Nm\n', max_efforts(4));

% --- Plot on a single figure with separate subplots ---
figure(3); clf;
set(gcf, 'Color', 'w');

% Motor 1: Base Rotation (Torque)
subplot(4, 1, 1);
plot(T, Fq_history(1,:), 'r', 'LineWidth', 1.5);
grid on;
ylabel('\tau_1 [Nm]');
title('Dynamic Motor Efforts Along the Trajectory');

% Motor 2: Prismatic Translation (Force)
subplot(4, 1, 2);
plot(T, Fq_history(2,:), 'b', 'LineWidth', 1.5);
grid on;
ylabel('F_2 [N]');

% Motor 3: Elbow Rotation (Torque)
subplot(4, 1, 3);
plot(T, Fq_history(3,:), 'g', 'LineWidth', 1.5);
grid on;
ylabel('\tau_3 [Nm]');

% Motor 4: Wrist Rotation (Torque)
subplot(4, 1, 4);
plot(T, Fq_history(4,:), 'y', 'LineWidth', 1.5);
grid on;
ylabel('\tau_4 [Nm]');
xlabel('Time [s]');