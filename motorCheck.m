% motorCheck.m - motor sizing: run AFTER projectA_main (uses T, Q, Qp, Qpp, atip, L_dim, alpha)
% Project sheet: "compute the motor torque of all the motors and check it with respect to
% the motor performance curve" + "performance at varying payload (0.2 kg to 1 kg) and cycle time"
%   - joint efforts for 2 materials x 2 payloads (your InverseDynamics + payload at the TCP)
%   - efforts reduced to the motor shaft, checked against the HG-KR23 data
%   - minimum reduction ratio needed, and effect of a shorter cycle time
% [ASSUMPTION] = value to replace with your design / catalogue data.
if ~exist('Qpp', 'var') || ~exist('atip', 'var'), error('Run projectA_main first.'); end

%% 1. DATA ------------------------------------------------------------------
matName = {'ABS'};
rhoMat  = 1050;                 % typical densities [kg/m^3]
% volumes of the 4 moving bodies [m^3] (SolidWorks: Evaluate > Mass Properties)
%   body 1 (q1): Part2 + 2 x Part3 | body 2 (q2): Part4 | body 3 (q3): Part5 | body 4 (q4): Part6
Vb = [12.915 + 2*2.174; 11.637; 19.538; 12.596]/2700;   % [ASSUMPTION] current SOLID parts -> update after lightening
payloads = [0.2 1.0];                   % [kg] project sheet
kT       = [1 0.5 0.25];                % cycle time x1, x0.5, x0.25 (same path, faster)

% motor: HG-KR23 (200 W) on every joint -> change per joint if you pick a smaller one
Tn   = 0.64*ones(4,1);                  % rated torque   [N m]
Tmax = 1.9 *ones(4,1);                  % maximum torque [N m]
nN   = 3000;  nMax = 4500;              % rated / maximum speed [rpm]
Jm   = 0.375e-4*ones(4,1);              % [ASSUMPTION] rotor inertia [kg m^2] -> check the catalogue
eta  = 0.9*ones(4,1);                   % [ASSUMPTION] efficiency of your own reduction stages
lead = 0.005;                           % [ASSUMPTION] lead of the q2 screw [m/rev]
iRed = [30; 2*pi/lead; 100; 20];        % [ASSUMPTION] reduction ratios (q2 in rad of motor per m)
jName = {'q1 (base)', 'q2 (prismatic)', 'q3 (elbow)', 'q4 (wrist)'};

%% 2. SIGN CONVENTION OF InverseDynamics (gravity must push q2 the same way) -
F0  = InverseDynamics(Q(:,1), zeros(4,1), zeros(4,1), L_dim, alpha, massesFor(rhoMat(1), Vb));
sgn = sign(F0(2));  if sgn == 0, sgn = 1; end       % payload term J'*F written with the same sign

%% 3. MATERIALS x PAYLOADS (planned cycle time) ---------------------------
nCase = numel(matName)*numel(payloads);  R = struct([]);  c = 0;
for im = 1:numel(matName)
    for ip = 1:numel(payloads)
        c = c + 1;
        Mass = massesFor(rhoMat(im), Vb);
        [Fq, Qps, Qpps] = jointEfforts(Q, Qp, Qpp, atip, L_dim, alpha, Mass, payloads(ip), 1, sgn);
        r = struct('name', sprintf('%s, %.1f kg', matName{im}, payloads(ip)), ...
                   'Fq', Fq, 'mass', [Mass.m1 Mass.m2 Mass.m3 Mass.m4]);
        r = motorSide(r, Fq, Qps, Qpps, T, iRed, eta, Jm);
        if c == 1, R = r; else, R(c) = r; end
    end
end

fprintf('\n=== JOINT SIDE: peak / RMS effort and minimum reduction ratio (HG-KR23) ===\n');
fprintf('(q1, q3, q4 in N m; q2 in N; for q2 the "ratio" is shown as the maximum screw lead)\n');
for c = 1:nCase
    fprintf('\n%s   (moving masses: %s kg)\n', R(c).name, num2str(R(c).mass, '%6.2f'));
    for j = 1:4
        Fpk  = max(abs(R(c).Fq(j,:)));
        Frms = rmsT(T, R(c).Fq(j,:));
        iMin = max(Fpk/(eta(j)*Tmax(j)), Frms/(eta(j)*Tn(j)));       % peak and continuous
        iSpd = (nMax*2*pi/60)/max(abs(Qp(j,:)) + eps);               % speed limit
        if j == 2
            fprintf('  %-15s peak %8.2f  RMS %8.2f  -> lead <= %6.1f mm (speed: lead >= %6.3f mm)\n', ...
                    jName{j}, Fpk, Frms, 1e3*2*pi/iMin, 1e3*2*pi/iSpd);
        else
            fprintf('  %-15s peak %8.2f  RMS %8.2f  -> ratio >= %6.1f   (speed: ratio <= %8.0f)\n', ...
                    jName{j}, Fpk, Frms, iMin, iSpd);
        end
    end
end

fprintf('\n=== MOTOR SIDE with your ratios iRed = [%s] (q2 lead %.1f mm) ===\n', ...
        num2str(iRed([1 3 4]).', '%g '), 1e3*lead);
fprintf('%-18s %-15s %9s %9s %9s  %s\n', 'case', 'joint', 'Tpk[Nm]', 'Trms[Nm]', 'n[rpm]', 'check');
for c = 1:nCase
    for j = 1:4
        ok = R(c).Tpk(j) <= Tmax(j) && R(c).Trms(j) <= Tn(j) && R(c).npk(j) <= nMax;
        fprintf('%-18s %-15s %9.3f %9.3f %9.0f  %s\n', R(c).name, jName{j}, ...
                R(c).Tpk(j), R(c).Trms(j), R(c).npk(j), passTxt(ok));
    end
end

%% 4. CYCLE TIME (payload 1 kg) ------------------------------------------------
TpkT = zeros(numel(matName), numel(kT), 4);  npkT = TpkT;
for im = 1:numel(matName)
    Mass = massesFor(rhoMat(im), Vb);
    for is = 1:numel(kT)
        [Fq, Qps, Qpps] = jointEfforts(Q, Qp, Qpp, atip, L_dim, alpha, Mass, max(payloads), kT(is), sgn);
        r = motorSide(struct(), Fq, Qps, Qpps, T*kT(is), iRed, eta, Jm);
        TpkT(im,is,:) = r.Tpk;  npkT(im,is,:) = r.npk;
    end
end
fprintf('\n=== CYCLE TIME (payload %.1f kg): peak motor torque [N m] / peak speed [rpm] ===\n', max(payloads));
for im = 1:numel(matName)
    for is = 1:numel(kT)
        fprintf('%-10s cycle x%-5.2f (%5.1f s):', matName{im}, kT(is), kT(is)*T(end));
        fprintf('  %6.3f/%5.0f', [squeeze(TpkT(im,is,:)).'; squeeze(npkT(im,is,:)).']);
        fprintf('\n');
    end
end

%% 5. FIGURES ------------------------------------------------------------------
col = lines(nCase);
lab = {'\tau_1 [N m]', 'F_2 [N]', '\tau_3 [N m]', '\tau_4 [N m]'};
figure(20); clf; set(gcf, 'Color', 'w', 'Name', 'Joint efforts');
for j = 1:4
    subplot(4,1,j); hold on; grid on;
    for c = 1:nCase, plot(T, R(c).Fq(j,:), 'Color', col(c,:), 'LineWidth', 1.2); end
    ylabel(lab{j});
    if j == 1, title('Joint efforts: materials x payloads'); legend({R.name}, 'Location', 'bestoutside'); end
end
xlabel('t [s]');

figure(21); clf; set(gcf, 'Color', 'w', 'Name', 'Motor torque-speed check');
for j = 1:4
    subplot(2,2,j); hold on; grid on;
    hE = [plot([0 nN nN], [Tn(j) Tn(j) 0], 'g-', 'LineWidth', 1.5), ...        % continuous zone (simplified)
          plot([0 nMax nMax], [Tmax(j) Tmax(j) 0], 'r--', 'LineWidth', 1.5)]; % short-duration zone (simplified)
    hS = gobjects(1, nCase);
    for c = 1:nCase
        plot(abs(R(c).n(j,:)), abs(R(c).Tm(j,:)), '.', 'Color', col(c,:), 'MarkerSize', 4);
        hS(c) = plot(R(c).nrms(j), R(c).Trms(j), 'p', 'Color', col(c,:), 'MarkerFaceColor', col(c,:), 'MarkerSize', 10);
    end
    if j == 1, legend([hE hS], [{'continuous (rated)', 'maximum'}, {R.name}], 'Location', 'best', 'FontSize', 7); end
    xlabel('motor speed [rpm]'); ylabel('motor torque [N m]'); title(jName{j});
    xlim([0 1.1*nMax]); ylim([0 1.2*max(Tmax(j), max(arrayfun(@(r) r.Tpk(j), R)))]);
end
sgtitle('HG-KR23: operating points (dots) and RMS points (stars)');

figure(22); clf; set(gcf, 'Color', 'w', 'Name', 'Cycle time');
for j = 1:4
    subplot(1,4,j); hold on; grid on;
    for im = 1:numel(matName), plot(kT*T(end), squeeze(TpkT(im,:,j)), 'o-', 'LineWidth', 1.2); end
    yline(Tmax(j), 'r--');  xlabel('cycle time [s]'); ylabel('peak motor torque [N m]'); title(jName{j});
end
legend([matName, {'T_{max}'}], 'Location', 'best');

%% LOCAL FUNCTIONS ---------------------------------------------------------------
function Mass = massesFor(rho, Vb)
% masses of the 4 moving bodies for one material
Mass.m1 = rho*Vb(1);  Mass.m2 = rho*Vb(2);  Mass.m3 = rho*Vb(3);  Mass.m4 = rho*Vb(4);
end

function [Fq, Qps, Qpps] = jointEfforts(Q, Qp, Qpp, atip, L_dim, alpha, Mass, mp, s, sgn)
% joint efforts along the path run with cycle time x s:
% robot (your InverseDynamics) + payload mp at the TCP (J' * force on the payload)
Qps = Qp/s;  Qpps = Qpp/s^2;  a = atip/s^2;               % same path, time scaled
M = size(Q, 2);  Fq = zeros(4, M);
for k = 1:M
    J   = Robot4DOF_jac(Q(:,k), L_dim, alpha);
    Fr  = InverseDynamics(Q(:,k), Qps(:,k), Qpps(:,k), L_dim, alpha, Mass);
    Fp  = mp*(a(:,k) + [0; 0; 9.81]);                     % force needed on the payload (gravity + inertia)
    Fq(:,k) = Fr(:) + sgn*J(1:3,:).'*Fp;
end
end

function r = motorSide(r, Fq, Qps, Qpps, t, iRed, eta, Jm)
% efforts at the motor shaft: load / (i*eta) + rotor inertia
r.Tm   = Fq./(iRed.*eta) + Jm.*iRed.*Qpps;                % [N m]
r.n    = iRed.*Qps*60/(2*pi);                             % [rpm]
r.Tpk  = max(abs(r.Tm), [], 2);
r.npk  = max(abs(r.n), [], 2);
r.Trms = zeros(4,1);  r.nrms = zeros(4,1);
for j = 1:4, r.Trms(j) = rmsT(t, r.Tm(j,:));  r.nrms(j) = rmsT(t, r.n(j,:)); end
end

function v = rmsT(t, x)
% RMS over the cycle (time weighted)
v = sqrt(trapz(t, x.^2)/(t(end) - t(1)));
end

function s = passTxt(ok)
if ok, s = 'OK'; else, s = 'NOT OK'; end
end