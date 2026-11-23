% projectA_mapping.m - Project A: the two "mapping" requests of the project sheet
%   (1) Map the effective workspace (with the mechanical limits of the design)
%   (2) Map the "quality" of the kinematics in the workspace
% Index used (course "Manipulability ellipses", Giberti):
%   M = J*J'  ->  sigma_i = sqrt(eig(M))                       (slides 13, 22)
%   k_inv = sigma_min/sigma_max  (inverse conditioning number, 0..1,
%           =1 isotropy, =0 singularity, > 0.2 acceptable)      (slide 23)
%   w     = sqrt(det(J*J')) = prod(sigma_i)  (ellipsoid volume, ~ det J)  (slide 22)
%   Velocity manipulability ellipses drawn on the letter                   (slides 10-14)
% Uses YOUR functions: Robot4DOF_inv, Robot4DOF_jac (same calls as projectA_main.m).
% Units: m, rad.   [ASSUMPTION] = value you must replace with data of your CAD.
clear; clc; close all;

%% 0. DATA - keep identical to projectA_main.m ----------------------------
R0 = 0.85;  Amp = 0.034;  Hs = 0.85;  kr = 2.5;
alpha = pi/6;
L1 = 0.29162 - 0.05103*sin(alpha);  L3 = 0.306;  L4 = 0.11928 + 0.10853/pi;
L_dim = [L1 L3 L4];
zc = 0.50;  zb = zc - Hs/2;  th0 = 0;
rS   = @(z) R0 + Amp*sin(2*pi*kr*(z - zb)/Hs);
drS  = @(z) Amp*2*pi*kr/Hs*cos(2*pi*kr*(z - zb)/Hs);
phiS = @(z) -atan(drS(z));                                % tool pitch = surface normal

%% 1. MECHANICAL LIMITS OF THE DESIGN -------------------------------------
q2min = 0.28 + 0.05103 + 0.2275 - 0.14;                   % prismatic folded   [m]
q2max = 0.28 + 0.05103 + 0.2275;                          % prismatic extended [m]
q1Lim = [-pi  pi];                                        % [ASSUMPTION] base rotation range [rad]
q3Lim = [-120 120]*pi/180;                              % [ASSUMPTION] elbow range, from the CAD [rad]
q4Lim = [-120 120]*pi/180;                              % [ASSUMPTION] wrist range, from the CAD [rad]
lim = [q1Lim; q2min q2max; q3Lim; q4Lim];                 % 4 x 2

%% 2. INDEX SETTINGS ------------------------------------------------------
rowsJ = 1:3;      % rows of J used: x,y,z (the ones verified in main, line 118).
                  % 1:4 would add the pitch row: rad vs m -> normalise first (slide 24)
Tn    = eye(4);   % normalisation of the joints, J*Tn (slides 24-26); eye = none
kMin  = 0.2;      % acceptable threshold of the inverse conditioning number (slide 23)

%% 3. TASK-SPACE GRID in the half-plane (rho, z): q1 only rotates the picture
dg   = 0.01;
rhoV = 0.20:dg:1.00;   zV = -0.10:dg:1.10;
[RHO, Z] = meshgrid(rhoV, zV);                            % rows = z, columns = rho
nG = numel(RHO);

%% 4. MAP 1 - EFFECTIVE WORKSPACE ------------------------------------------
% (a) any tool pitch: point reachable if at least one pitch phi works
phiSet = (-90:15:90)*pi/180;
WS_any = false(size(RHO));
for i = 1:nG
    for p = 1:numel(phiSet)
        [~, ok] = ikLim([RHO(i); 0; Z(i); phiSet(p)], L_dim, alpha, lim);
        if ok, WS_any(i) = true; break; end
    end
end
% (b) task orientation: tool normal to the surface, phi = phiS(z)
WS_task = false(size(RHO));  K = nan(size(RHO));  W = nan(size(RHO));
for i = 1:nG
    [Q, ok] = ikLim([RHO(i); 0; Z(i); phiS(Z(i))], L_dim, alpha, lim);
    if ~ok, continue; end
    WS_task(i) = true;
    [K(i), W(i)] = qualityIdx(Q, L_dim, alpha, rowsJ, Tn);
end

% wall of the surface: is every point of r(z) inside the task workspace?
zw = linspace(zb, zb+Hs, 400);
okW = false(size(zw));
for i = 1:numel(zw)
    [~, okW(i)] = ikLim([rS(zw(i)); 0; zw(i); phiS(zw(i))], L_dim, alpha, lim);
end
fprintf('Workspace area (rho-z half plane): any pitch = %.3f m^2, tool normal to surface = %.3f m^2\n', ...
        nnz(WS_any)*dg^2, nnz(WS_task)*dg^2);
fprintf('Surface wall reachable with tool normal: %.1f %% of its height (z = %.3f .. %.3f m)\n', ...
        100*mean(okW), min(zw(okW)), max(zw(okW)));

%% 5. THE LETTER (same nodes as main, section 2) ---------------------------
V = [-30  -18.75  18.75    0     30;
     -40  -10    -10      40    -40]*1e-3;
order = [1 2 3 4 2 3 5];
dsNode = 0.002;
x = []; y = [];
for i = 1:numel(order)-1
    a = V(:,order(i)); b = V(:,order(i+1));
    n = max(1, round(norm(b - a)/dsNode));
    s = (0:n-1)/n;
    x = [x, a(1) + (b(1)-a(1))*s];  %#ok<AGROW>
    y = [y, a(2) + (b(2)-a(2))*s];  %#ok<AGROW>
end
x(end+1) = V(1,order(end));  y(end+1) = V(2,order(end));
z  = zc + y;  th = th0 + x/R0;  r = rS(z);
S_let = [r.*cos(th); r.*sin(th); z; phiS(z)];
N = size(S_let, 2);

okL = false(1,N);  KL = nan(1,N);  WL = nan(1,N);  QL = nan(4,N);
for i = 1:N
    [Q, ok] = ikLim(S_let(:,i), L_dim, alpha, lim);
    okL(i) = ok;
    if ok, QL(:,i) = Q;  [KL(i), WL(i)] = qualityIdx(Q, L_dim, alpha, rowsJ, Tn); end
end
fprintf('Letter: %d/%d nodes inside the effective workspace\n', nnz(okL), N);
fprintf('Letter: k_inv min = %.3f, mean = %.3f, max = %.3f (acceptable > %.1f)\n', ...
        min(KL), mean(KL,'omitnan'), max(KL), kMin);

%% 6. FIGURE 1 - effective workspace --------------------------------------
figure(1); clf; set(gcf,'Color','w');
subplot(1,2,1); hold on; grid on; axis equal;
layer = double(WS_any) + double(WS_task);                 % 0 outside, 1 any pitch, 2 + normal to surface
imagesc(rhoV, zV, layer);  set(gca,'YDir','normal');
colormap([1 1 1; 0.85 0.85 0.85; 0.55 0.75 0.95]);  caxis([0 2]);
zs = linspace(zb, zb+Hs, 400);
hS = plot(rS(zs), zs, 'r', 'LineWidth', 1.5);
hL = plot(sqrt(S_let(1,:).^2 + S_let(2,:).^2), S_let(3,:), 'k.', 'MarkerSize', 4);
plot([0 0], [0 L1], 'k', 'LineWidth', 3);
xlim([0 1]); ylim([-0.1 1.1]);
xlabel('\rho [m]'); ylabel('z [m]');
title('Effective workspace (half plane, q_1 rotates it)');
legend([hS hL], {'surface r(z)', 'letter A'}, 'Location', 'southeast');
text(0.22, 1.05, 'grey: any pitch   blue: tool normal to surface', 'FontSize', 8);

subplot(1,2,2); hold on; grid on; axis equal; view(35, 20);
[iz, ir] = find(WS_task);
thv = linspace(q1Lim(1), q1Lim(2), 40);
for t = 1:numel(thv)
    scatter3(rhoV(ir)*cos(thv(t)), rhoV(ir)*sin(thv(t)), zV(iz), 2, [0.35 0.6 0.9], '.');
end
[Zs, Ts] = meshgrid(linspace(zb, zb+Hs, 60), linspace(-0.6, 0.6, 60));
surf(rS(Zs).*cos(Ts), rS(Zs).*sin(Ts), Zs, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'r');
plot3(S_let(1,:), S_let(2,:), S_let(3,:), 'k.', 'MarkerSize', 4);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Workspace with tool normal to surface (3D)');

%% 7. FIGURE 2 - quality maps ----------------------------------------------
figure(2); clf; set(gcf,'Color','w');
subplot(1,2,1); hold on; grid on; axis equal;
pcolor(RHO, Z, K); shading flat;  colorbar;  caxis([0 0.5]);
contour(RHO, Z, K, [kMin kMin], 'k', 'LineWidth', 1.5);
plot(rS(zs), zs, 'r', 'LineWidth', 1.5);
plot(sqrt(S_let(1,:).^2 + S_let(2,:).^2), S_let(3,:), 'k.', 'MarkerSize', 4);
xlabel('\rho [m]'); ylabel('z [m]'); xlim([0 1]); ylim([-0.1 1.1]);
title('k_{inv} = \sigma_{min}/\sigma_{max}  (black line: 0.2)');

subplot(1,2,2); hold on; grid on; axis equal;
pcolor(RHO, Z, W); shading flat;  colorbar;
plot(rS(zs), zs, 'r', 'LineWidth', 1.5);
plot(sqrt(S_let(1,:).^2 + S_let(2,:).^2), S_let(3,:), 'k.', 'MarkerSize', 4);
xlabel('\rho [m]'); ylabel('z [m]'); xlim([0 1]); ylim([-0.1 1.1]);
title('w = sqrt(det(JJ^T))  (ellipsoid volume)');

%% 8. FIGURE 3 - quality along the letter + manipulability ellipses --------
figure(3); clf; set(gcf,'Color','w');
sArc = [0, cumsum(sqrt(sum(diff(S_let(1:3,:),1,2).^2,1)))];
subplot(1,2,1); hold on; grid on;
plot(sArc, KL, 'b', 'LineWidth', 1.5);  yline(kMin, 'r--', '0.2');
xlabel('path along the letter [m]'); ylabel('k_{inv}'); ylim([0 1]);
title('Kinematic quality along the letter');

subplot(1,2,2); hold on; grid on; axis equal;
plot(x, y + zc, 'k-', 'LineWidth', 1);
iE   = round(linspace(1, N, 12));
scl  = 0.012;                                              % display scale of the ellipses [1/(m/s per rad/s)]
tt   = linspace(0, 2*pi, 60);
for i = iE
    if ~okL(i), continue; end
    J    = Robot4DOF_jac(QL(:,i), L_dim, alpha);
    eTh  = [-sin(th(i)); cos(th(i)); 0];                  % tangential direction on the surface
    Jloc = [eTh.'; 0 0 1] * J(1:3,:) * Tn;                % (tangential, vertical) x joints
    [Ve, De] = eig(Jloc*Jloc.');                          % M = J*J'  (slide 13)
    e = Ve*sqrt(max(De,0))*[cos(tt); sin(tt)];            % image of the unit joint-velocity circle
    plot(x(i) + scl*e(1,:), z(i) + scl*e(2,:), 'b');
    plot(x(i), z(i), 'r.', 'MarkerSize', 8);
end
xlabel('tangential [m]'); ylabel('z [m]');
title('Velocity manipulability ellipses on the letter (scaled)');

%% LOCAL FUNCTIONS ---------------------------------------------------------
function [Q, ok] = ikLim(S, L_dim, alpha, lim)
% IK of the project + feasibility: real solution inside the joint limits
Q  = Robot4DOF_inv(S, L_dim, alpha);
ok = all(abs(imag(Q)) < 1e-12) && all(isfinite(Q));
Q  = real(Q);
ok = ok && all(Q >= lim(:,1) - 1e-9) && all(Q <= lim(:,2) + 1e-9);
end

function [kinv, w] = qualityIdx(Q, L_dim, alpha, rowsJ, Tn)
% inverse conditioning number and ellipsoid volume of the Jacobian
J  = Robot4DOF_jac(Q, L_dim, alpha);
sv = svd(J(rowsJ,:) * Tn);                                % sigma_i = sqrt(eig(J*J'))
kinv = sv(end)/sv(1);
w    = prod(sv);
end