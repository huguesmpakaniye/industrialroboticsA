% projectA_main.m - Project A: letter "A" on the sinusoidal cylinder
% Joint trajectories as required by the project sheet (p.4):
%   (1) home -> first point of the A : CIRCULAR path,  cubic splines in joint space
%   (2) letter A                     : lines and parabolas in joint space
%   (3) last point of the A -> home  : LINEAR path,    cubic splines in joint space
% Uses: Robot4DOF_inv, ForwardKinematics, Robot4DOF_jac, InverseDynamics,
%       PlotRobot4DOF, LinesParabolas, Circle3pt, liveTrace (Simulink),
%       motorSizing, performanceAnalysis (section 11)
% Simscape model: simulink_multibody_SCARA.slx (CAD parts in CAD/, density rhoABS)
% Separate script: projectA_mapping.m (workspace + kinematic quality maps)
% Units: m, rad, s.   All numbers marked [ASSUMPTION] are yours to change.
clear; clc; close all;
cd(fileparts(mfilename('fullpath')));   % work from the repo folder (CAD/ paths of the Simscape model)

%% 0. SURFACE AND ROBOT --------------------------------------------------
R0 = 0.72;  Amp = 0.04*R0;  Hs = R0;  kr = 2.5;          % A/R0 = 0.04, H/R0 = 1.0
alpha = pi/6;                                            % prismatic axis inclination (checked with Simscape)
L1 = 0.320 + 0.055;                                      % part1 height + Part2 trunnion pins
L3 = 0.245;                                              % Part5: 259.4 - 14.4 mm
L4 = 0.180;                                              % Part6: wrist axis -> end face
L_dim = [L1 L3 L4];
q2min = 0.1125 + 0.1775;  q2max = 0.2525 + 0.1775;       % 0.290 .. 0.430 m
zc = 0.60;  zb = zc - Hs/2;  th0 = 0;                    % letter centre 0.60 m, surface 0.24 -> 0.96 m
rS   = @(z) R0 + Amp*sin(2*pi*kr*(z - zb)/Hs);            % surface radius r(z)
drS  = @(z) Amp*2*pi*kr/Hs*cos(2*pi*kr*(z - zb)/Hs);      % r'(z)
phiS = @(z) -atan(drS(z));                                % tool pitch = surface normal
S_home = [rS(zc) - 0.050; 0; zc; 0];                      % [ASSUMPTION] home: 50 mm inside the wall, tool horizontal
Q_home = Robot4DOF_inv(S_home, L_dim, alpha);

%% 1. TIMING / LIMITS ----------------------------------------------------
vTool  = 0.010;                 % tool speed on the letter [m/s]
dsNode = 0.002;                 % node spacing along the letter [m]
accMax = [0.5; 1; 10; 20];      % blend accelerations [rad/s2; m/s2; rad/s2; rad/s2] (q1 moves only +-2 deg)
Tcirc  = 2.0;  nCirc = 13;      % stage 1: duration, spline nodes on the circle
Tret   = 3.0;  nLin  = 9;       % stage 3: duration, spline nodes on the straight line
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
assert(all(Q_let(2,:) >= q2min & Q_let(2,:) <= q2max), ...
       'q2 out of stroke on the letter (%.3f..%.3f)', min(Q_let(2,:)), max(Q_let(2,:)));

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

%% 4. STAGE 1 - CIRCULAR APPROACH: cubic splines in joint space ------------
P1 = S_home(1:3);  P3 = S_let(1:3,1);                     % home -> first point of the A
Mid = (P1 + P3)/2;  ch = P3 - P1;
up = [0;0;1];  perp = up - (up.'*ch)/(ch.'*ch)*ch;  perp = perp/norm(perp);
P2 = Mid + 0.010*perp;                                    % [ASSUMPTION] 10 mm bulge keeps the tip off the wall
Parc = Circle3pt(P1, P2, P3, nCirc);
S_arc = [Parc; linspace(S_home(4), S_let(4,1), nCirc)];
Q_arc = zeros(4, nCirc);
for i = 1:nCirc, Q_arc(:,i) = Robot4DOF_inv(S_arc(:,i), L_dim, alpha); end
assert(isreal(Q_arc), 'Circular approach not reachable');
tc  = linspace(0, Tcirc, nCirc);
[QC, QpC, QppC, tC] = jointSpline(tc, Q_arc, dt);         % clamped: zero velocity at both ends

%% 5. STAGE 3 - LINEAR RETURN: cubic splines in joint space ----------------
% nodes taken on the straight line (task space) last point of the A -> home,
% then the joints are interpolated with clamped cubic splines (as stage 1)
S_end = S_let(:,end);
sLin  = linspace(0, 1, nLin);
S_lin = S_end*(1 - sLin) + S_home*sLin;                   % 4 x nLin: position and pitch on the line
Q_lin = zeros(4, nLin);
for i = 1:nLin, Q_lin(:,i) = Robot4DOF_inv(S_lin(:,i), L_dim, alpha); end
assert(isreal(Q_lin), 'Linear return not reachable');
tl  = linspace(0, Tret, nLin);
[QR, QpR, QppR, tR] = jointSpline(tl, Q_lin, dt);         % clamped: zero velocity at both ends

%% 6. ASSEMBLE -------------------------------------------------------------
T   = [tC, Tcirc + tL(2:end).', Tcirc + tL(end) + tR(2:end)];
Q   = [QC,   QL(:,2:end),   QR(:,2:end)];
Qp  = [QpC,  QpL(:,2:end),  QpR(:,2:end)];
Qpp = [QppC, QppL(:,2:end), QppR(:,2:end)];
assert(all(Q(2,:) >= q2min & Q(2,:) <= q2max), ...
       'q2 out of stroke (%.3f..%.3f)', min(Q(2,:)), max(Q(2,:)));
iArc = T <= Tcirc;                                        % stage 1 samples
iLet = (T >= Tcirc) & (T <= Tcirc + tL(end));             % stage 2 samples
iRet = T >= Tcirc + tL(end);                              % stage 3 samples

%% 7. VERIFY WITH FORWARD KINEMATICS --------------------------------------
M = numel(T);  tip = zeros(3,M);  vtip = zeros(3,M);  atip = zeros(3,M);
for k = 1:M
    [M04, W04, H04] = ForwardKinematics(Q(:,k), Qp(:,k), Qpp(:,k), L_dim, alpha);
    p = M04(:,4);  tip(:,k) = p(1:3);
    vv = W04*p;  vtip(:,k) = vv(1:3);
    aa = H04*p;  atip(:,k) = aa(1:3);
end
kmid = round(M/2);
J = Robot4DOF_jac(Q(:,kmid), L_dim, alpha);
[~, Wm] = ForwardKinematics(Q(:,kmid), Qp(:,kmid), Qpp(:,kmid), L_dim, alpha);
vm = Wm*[tip(:,kmid); 1];
fprintf('Velocity check |W*p - J*qp| = %.2e (must be ~0)\n', norm(vm(1:3) - J(1:3,:)*Qp(:,kmid)));

rho  = sqrt(tip(1,:).^2 + tip(2,:).^2);
gap  = rS(tip(3,:)) - rho;                                % >0: tip is inside the wall, 0: on the wall
fprintf('On the letter: max |radial error| = %.2f mm\n', 1e3*max(abs(gap(iLet))));
fprintf('Min wall gap on the arc (excluding last 0.2 s) = %.1f mm\n', ...
        1e3*min(gap(T>0.05 & T<Tcirc-0.2)));
fprintf('Min wall gap on the return (excluding first 0.2 s) = %.1f mm\n', ...
        1e3*min(gap(T > Tcirc + tL(end) + 0.2)));
% stage 1: distance of the tip from the circle  /  stage 3: from the straight line
Cc   = circleCentre(P1, P2, P3);  Rarc = norm(P1 - Cc);
dArc = abs(vecnorm(tip(:,iArc) - Cc, 2, 1) - Rarc);
uL   = (S_home(1:3) - S_end(1:3))/norm(S_home(1:3) - S_end(1:3));
dP   = tip(:,iRet) - S_end(1:3);
dLin = vecnorm(dP - uL*(uL.'*dP), 2, 1);
fprintf('Arc: max distance from the circle = %.3f mm   Return: max distance from the line = %.3f mm\n', ...
        1e3*max(dArc), 1e3*max(dLin));
fprintf('Peak joint speeds   : %s\n', num2str(max(abs(Qp),[],2).', '%8.3f'));
fprintf('Peak joint accel.   : %s\n', num2str(max(abs(Qpp),[],2).', '%8.3f'));
fprintf('Joint ranges        : q1 [%.2f %.2f]  q2 [%.3f %.3f]  q3 [%.2f %.2f]  q4 [%.2f %.2f]\n', ...
        min(Q(1,:)),max(Q(1,:)),min(Q(2,:)),max(Q(2,:)),min(Q(3,:)),max(Q(3,:)),min(Q(4,:)),max(Q(4,:)));

%% 8. PLOTS ----------------------------------------------------------------
figure(1); clf; hold on; grid on; axis equal;
[Zs, Ts] = meshgrid(linspace(zb, zb+Hs, 60), linspace(-0.6, 0.6, 60));
surf(rS(Zs).*cos(Ts), rS(Zs).*sin(Ts), Zs, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot3(tip(1,iArc), tip(2,iArc), tip(3,iArc), 'b', 'LineWidth', 1.5);
plot3(tip(1,iLet), tip(2,iLet), tip(3,iLet), 'r', 'LineWidth', 1.5);
plot3(tip(1,iRet), tip(2,iRet), tip(3,iRet), 'g', 'LineWidth', 1.5);
plot3(S_let(1,:), S_let(2,:), S_let(3,:), 'k.', 'MarkerSize', 4);
PlotRobot4DOF(Q_home, L_dim, alpha, 1, 'b', '-', 'Project A trajectory');
PlotRobot4DOF(Q_let(:,1), L_dim, alpha, 1, 'g', '-', 'Project A trajectory');
xlim([-0.1 1.0]); ylim([-0.5 0.5]); zlim([0 1]);
figure(2); clf;
sp = {'q_1 [rad]','q_2 [m]','q_3 [rad]','q_4 [rad]'};
for j = 1:4
    subplot(4,3,3*j-2); plot(T, Q(j,:));   grid on; ylabel(sp{j}); if j==1, title('position'); end
    subplot(4,3,3*j-1); plot(T, Qp(j,:));  grid on;               if j==1, title('velocity');  end
    subplot(4,3,3*j);   plot(T, Qpp(j,:)); grid on;               if j==1, title('acceleration'); end
end
xlabel('t [s]');

%% 9. DYNAMIC ANALYSIS (Inverse Dynamics) ----------------------------------
% Material: ABS (project sheet: aluminium alloy or ABS).
% rhoABS is also used as Density in the File Solid blocks of the Simscape model,
% so Matlab and Simscape have the SAME masses (needed for the torque comparison).
rhoABS = 1050;                                            % [kg/m^3]
Vb = [12.915 + 2*2.174; 11.637; 19.538; 12.596]/2700;     % volumes of the moving bodies from the CAD [m^3]
Masses.m1 = rhoABS*Vb(1);   % Part2 + 2 x Part3 (rotate with q1)   ~6.71 kg
Masses.m2 = rhoABS*Vb(2);   % Part4 (slider)                       ~4.53 kg
Masses.m3 = rhoABS*Vb(3);   % Part5 (link L3)                      ~7.60 kg
Masses.m4 = rhoABS*Vb(4);   % Part6 (end effector)                 ~4.90 kg
Masses.mot3 = 1.41;   % motor HG-KR23 + reducer of q3, on the elbow axis (Simscape point mass 1)
Masses.mot4 = 0.65;   % motor HG-KR053 + reducer of q4, at mid L3     (Simscape point mass 2)
% (InverseDynamics places every mass at its CAD centre of mass, as Simscape does)
fprintf('\nMoving masses (ABS): %s kg\n', num2str([Masses.m1 Masses.m2 Masses.m3 Masses.m4], '%6.2f'));

Fq_history = zeros(4, M);
for k = 1:M
    Fq_history(:, k) = InverseDynamics(Q(:,k), Qp(:,k), Qpp(:,k), L_dim, alpha, Masses);
end
max_efforts = max(abs(Fq_history), [], 2);
fprintf('\n--- DYNAMIC RESULTS (Motor Sizing) ---\n');
fprintf('Max Tau 1 (Base)        : %.2f Nm\n', max_efforts(1));
fprintf('Max Force 2 (Prismatic) : %.2f N\n',  max_efforts(2));
fprintf('Max Tau 3 (Elbow)       : %.2f Nm\n', max_efforts(3));
fprintf('Max Tau 4 (Wrist)       : %.2f Nm\n', max_efforts(4));

figure(3); clf; set(gcf, 'Color', 'w');
lab = {'\tau_1 [Nm]','F_2 [N]','\tau_3 [Nm]','\tau_4 [Nm]'};
for j = 1:4
    subplot(4,1,j); plot(T, Fq_history(j,:), 'LineWidth', 1.5); grid on; ylabel(lab{j});
    if j==1, title('Dynamic Motor Efforts Along the Trajectory'); end
end
xlabel('Time [s]');

%% 10. SIMSCAPE MULTIBODY: joint inputs = the planned trajectory -----------
% Each joint is driven in MOTION (Simulink-PS Converter with 1st and 2nd
% derivatives provided) by its position, velocity and acceleration.
fprintf('\n--- PREPARING SIMSCAPE SIMULATION ---\n');
% Simscape joint value = Q_sgn .* (Q - Q_off)
%   Q_off : value of the IK joint variable when the Simscape joint is at 0.
%           q2 (IK) = TOTAL length from the axis intersection to the elbow (0.29 .. 0.43 m),
%           p  (Simscape) = slider DISPLACEMENT from its assembled pose (0 .. 0.14 m)
%   Q_sgn : -1 if a Simscape joint moves the wrong way
Q_off = [pi/2; q2min; 0; 0];   % q1: the arm faces +y in Simscape when its joint is at 0 (checked)
Q_sgn = [1; 1; -1; 1];         % q3 turns the other way in Simscape (checked with checkJointZeros)
Q_sim = (Q - Q_off).*Q_sgn;
fprintf('Simscape prismatic input p = %.3f .. %.3f m  (rail: 0 .. %.2f m)\n', ...
        min(Q_sim(2,:)), max(Q_sim(2,:)), q2max - q2min);
assert(all(Q_sim(2,:) >= -1e-9 & Q_sim(2,:) <= q2max - q2min + 1e-9), ...
       'Prismatic input outside the rail: check Q_off(2) / Q_sgn(2)');
T_sim  = T(:);
Q_ts   = timeseries(Q_sim.', T_sim);                       % M x 4: q1 q2 q3 q4 as Simscape joints
Qp_ts  = timeseries((Qp.*Q_sgn).',  T_sim);
Qpp_ts = timeseries((Qpp.*Q_sgn).', T_sim);

% --- surface for Mechanics Explorer: cross-section of a "Revolved Solid" block (x-z plane, m) ---
eW = 0.005;  zS = linspace(zb, zb + Hs, 200).';           % 5 mm thick wall behind the drawing face r(z)
surfXS = [rS(flipud(zS)), flipud(zS);  rS(zS) + eW, zS];  % closed, counterclockwise

% --- live figure: liveTrace.m (called by the model) adds the TCP points ---
figure(10); clf; set(gcf, 'Color', 'w', 'Name', 'Live TCP trace');
subplot(1,2,1); hold on; grid on; axis equal; view(60,20);
surf(rS(Zs).*cos(Ts), rS(Zs).*sin(Ts), Zs, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot3(S_let(1,:), S_let(2,:), S_let(3,:), 'k:');
animatedline('Tag', 'TCPtrace3D', 'Color', 'r', 'LineWidth', 2);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); title('TCP on the surface (live)');
subplot(1,2,2); hold on; grid on; axis equal;
plot(R0*th, z, 'k:');                                     % planned letter, unrolled (R0*theta, z)
animatedline('Tag', 'TCPtrace2D', 'Color', 'r', 'LineWidth', 2, ...
             'UserData', struct('R0', R0, 'rS', rS, 'tol', 1e-3));   % drawn only when the pen touches the wall
xlim([-0.05 0.05]); ylim(zc + [-0.06 0.06]);
xlabel('R_0 \theta [m]'); ylabel('z [m]'); title('Letter drawn on the wall (live)');
drawnow;

mdl = 'simulink_multibody_SCARA';                         % model name WITHOUT extension
load_system(mdl);
try
    set_param(mdl, 'EnablePacing', 'on', 'PacingRate', 1);   % real time, to watch it live
catch
    warning('Set Run > Simulation Pacing by hand to watch the simulation in real time.');
end
disp('Starting 3D simulation...');
out = sim(mdl, 'StopTime', num2str(T(end)));
disp('Simulation finished!');

% --- check: Simscape TCP must follow the Matlab FK (otherwise: L, frames or Q_off) ---
[tS, tipS] = getSig(out, 'tipSim');
tipSi = interp1(tS, tipS, T_sim, 'linear', 'extrap');     % M x 3
errTCP = vecnorm(tipSi.' - tip, 2, 1);
[eMax, iMax] = max(errTCP);
fprintf('Simscape vs Matlab TCP: max error = %.2f mm, at t = %.2f s\n', 1e3*eMax, T(iMax));

figure(4); clf; hold on; grid on; axis equal;
surf(rS(Zs).*cos(Ts), rS(Zs).*sin(Ts), Zs, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot3(tipS(:,1), tipS(:,2), tipS(:,3), 'r', 'LineWidth', 1.5);
plot3(tip(1,:), tip(2,:), tip(3,:), 'b--');
legend('surface', 'Simscape', 'Matlab FK'); view(60,20);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

% --- torques Matlab vs Simscape (needs a To Workspace 'tauSim' = 4 actuator torques/force) ---
if any(strcmp(out.who, 'tauSim'))
    [tT, tauS] = getSig(out, 'tauSim');
    tauSi = interp1(tT, tauS, T_sim, 'linear', 'extrap');
    figure(5); clf; set(gcf, 'Color', 'w');
    for j = 1:4
        subplot(4,1,j); hold on; grid on;
        plot(T, Fq_history(j,:), 'b', 'LineWidth', 1.5);
        plot(T, Q_sgn(j)*tauSi(:,j), 'r--', 'LineWidth', 1.2);   % same sign convention as Matlab
        ylabel(lab{j}); if j==1, title('Motor efforts: Matlab vs Simscape'); legend('Matlab','Simscape'); end
    end
    xlabel('Time [s]');
end

%% 11. MOTOR SIZING AND PERFORMANCE ----------------------------------------
motorSizing;            % motors + reducers (course: dynamic sizing)
performanceAnalysis;    % payload 0.2-1 kg and cycle time

%% LOCAL FUNCTIONS ---------------------------------------------------------
function [Qs, Qps, Qpps, ts] = jointSpline(tn, Qn, dt)
% clamped cubic spline in joint space through the nodes Qn (4 x n) at times tn
pp   = spline(tn, [zeros(size(Qn,1),1), Qn, zeros(size(Qn,1),1)]);   % zero velocity at both ends
c    = pp.coefs;
ppd  = mkpp(pp.breaks, [3*c(:,1) 2*c(:,2) c(:,3)], pp.dim);
ppdd = mkpp(pp.breaks, [6*c(:,1) 2*c(:,2)], pp.dim);
ts   = 0:dt:tn(end);
if ts(end) < tn(end) - 1e-9, ts(end+1) = tn(end); end
Qs = ppval(pp, ts);  Qps = ppval(ppd, ts);  Qpps = ppval(ppdd, ts);
end

function C = circleCentre(P1, P2, P3)
% centre of the circle through three points
a = P1 - P3;  b = P2 - P3;  n = cross(a, b);
C = P3 + cross(dot(a,a)*b - dot(b,b)*a, n)/(2*dot(n,n));
end

function [t, d] = getSig(out, name)
% signal saved by a To Workspace block (Array or Timeseries format) -> t (N x 1), d (N x n)
s = out.get(name);
if isa(s, 'timeseries'), t = s.Time; d = squeeze(s.Data);
else,                    t = out.tout; d = squeeze(s); end
if size(d,1) ~= numel(t), d = d.'; end
[t, iu] = unique(t);  d = d(iu,:);                       % interp1 needs distinct times
end