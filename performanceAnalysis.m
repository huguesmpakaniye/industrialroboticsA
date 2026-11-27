% performanceAnalysis.m - project sheet: "Analyse the robot's performance for the task at
% varying payload (0.2 kg to 1 kg) and cycle time".
% Run AFTER projectA_main AND motorSizing (uses the chosen motors S(j).pick and ratios S(j).tau).
% Same path, run faster or slower: cycle time x s  ->  qdot/s, qddot/s^2, tip acceleration/s^2.
% Exact decomposition of the joint efforts (no approximation):
%   Cr*(mp, s) = G + D/s^2 + mp*(Pg + Pa/s^2)
%   G  = gravity of the robot,  D  = inertia + Coriolis of the robot at s = 1,
%   Pg = payload gravity per kg, Pa = payload inertia per kg at s = 1
% Motor side as in the course:  Cm = tau*Cr* + Jm*wd_r/tau ;  checks Cmax, Crms < Cn, w < w_max
if ~exist('S', 'var') || ~isfield(S, 'pick'), error('Run projectA_main then motorSizing first.'); end
if any([S.pick] == 0), error('motorSizing found no motor for one joint.'); end

%% 1. SETTINGS --------------------------------------------------------------------
mpList = 0.2:0.1:1.0;                    % payloads of the project sheet [kg]
sList  = logspace(0, log10(0.02), 400);  % cycle-time factors: 1 (planned) -> 0.02 (50x faster)
mpMaxScan = 0:0.1:10;                    % to find the largest payload at the planned cycle time
jName = {'q1 (base)', 'q2 (prismatic)', 'q3 (elbow)', 'q4 (wrist)'};
tc0 = T(end) - T(1);                     % planned cycle time [s]
k   = [S.pick];  tau = arrayfun(@(j) S(j).tau(k(j)), 1:4).';
tauOverride = [1/40; NaN; NaN; NaN];   % [ASSUMPTION] q1: 1/tau = 40 (one belt stage) for faster cycles
tau(~isnan(tauOverride)) = tauOverride(~isnan(tauOverride));
CnJ = Cn(k).';  CmxJ = CmMax(k).';  JmJ = Jm(k).';  wmxJ = wmMax(k).';

%% 2. THE FOUR TERMS (computed once) ------------------------------------------------
M  = numel(T);  G = zeros(4,M);  D = G;  Pg = G;  Pa = G;
for i = 1:M
    Fs = InverseDynamics(Q(:,i), zeros(4,1), zeros(4,1), L_dim, alpha, Masses);
    Fd = InverseDynamics(Q(:,i), Qp(:,i),    Qpp(:,i),   L_dim, alpha, Masses);
    J  = Robot4DOF_jac(Q(:,i), L_dim, alpha);
    G(:,i)  = Fs(:);   D(:,i) = Fd(:) - Fs(:);
    Pg(:,i) = J(1:3,:).'*[0; 0; 9.81];
    Pa(:,i) = J(1:3,:).'*atip(:,i);
end
sgn = sign(G(2,1));  if sgn == 0, sgn = 1; end       % same sign convention as InverseDynamics
Pg = sgn*Pg;  Pa = sgn*Pa;

%% 3. CHECKS FOR EVERY (payload, cycle time) ----------------------------------------
% utilisation u = max(Cpk/Cmax, Crms/Cn, w/w_max): the motor is OK while u < 1
nP = numel(mpList);  nS = numel(sList);
U  = zeros(4, nP, nS);  Upk = U;  Urms = U;  Uw = U;
for ip = 1:nP
    for is = 1:nS
        [U(:,ip,is), Upk(:,ip,is), Urms(:,ip,is), Uw(:,ip,is)] = util(mpList(ip), sList(is), ...
            G, D, Pg, Pa, Qp, Qpp, T, tau, CnJ, CmxJ, JmJ, wmxJ);
    end
end

% minimum cycle time: smallest s such that every faster-than-planned step down to it is OK
tcMin = nan(1, nP);  limJ = zeros(1, nP);  limC = strings(1, nP);
chk = ["peak torque", "RMS torque (thermal)", "speed"];
for ip = 1:nP
    okS = all(squeeze(U(:,ip,:)) < 1, 1);              % 1 x nS
    last = find(~okS, 1) - 1;                           % last OK index before the first failure
    if isempty(last), last = nS; end
    if last == 0, continue; end                         % already not OK at the planned cycle
    tcMin(ip) = sList(last)*tc0;
    if last < nS                                        % what fails just after
        [~, limJ(ip)] = max(U(:,ip,last+1));
        [~, c] = max([Upk(limJ(ip),ip,last+1), Urms(limJ(ip),ip,last+1), Uw(limJ(ip),ip,last+1)]);
        limC(ip) = chk(c);
    end
end

% largest payload at the planned cycle time
okMp = arrayfun(@(mp) all(util(mp, 1, G, D, Pg, Pa, Qp, Qpp, T, tau, CnJ, CmxJ, JmJ, wmxJ) < 1), mpMaxScan);
iBad = find(~okMp, 1);
if isempty(iBad), mpMax = mpMaxScan(end); else, mpMax = mpMaxScan(max(iBad-1, 1)); end
[uMax, jMax] = max(util(mpMax + 0.1, 1, G, D, Pg, Pa, Qp, Qpp, T, tau, CnJ, CmxJ, JmJ, wmxJ));

%% 4. PRINT ---------------------------------------------------------------------------
fprintf('\n=== PERFORMANCE: chosen motors and transmissions ===\n');
for j = 1:4
    if j == 2, r = sprintf('lead %.1f mm', 1e3*2*pi*tau(j)); else, r = sprintf('1/tau = %.1f', 1/tau(j)); end
    fprintf('%-15s %-9s %s\n', jName{j}, motName{k(j)}, r);
end
fprintf('\nPlanned cycle time: %.1f s\n', tc0);
fprintf('%-8s %-16s %-10s %s\n', 'payload', 'min cycle [s]', 'speed-up', 'first limit reached');
for ip = 1:nP
    if isnan(tcMin(ip)), fprintf('%5.1f kg  not OK at the planned cycle time\n', mpList(ip)); continue; end
    if limJ(ip) == 0, lim = 'none in the range studied'; else, lim = sprintf('%s, %s', jName{limJ(ip)}, limC(ip)); end
    fprintf('%5.1f kg  %10.1f       x%-8.1f %s\n', mpList(ip), tcMin(ip), tc0/tcMin(ip), lim);
end
fprintf('\nLargest payload at the planned cycle time: %.1f kg (then %s is the first to exceed, u = %.2f)\n', ...
        mpMax, jName{jMax}, uMax);
fprintf('Utilisation at planned cycle, 1 kg: %s  (q1..q4, must be < 1)\n', ...
        num2str(util(1.0, 1, G, D, Pg, Pa, Qp, Qpp, T, tau, CnJ, CmxJ, JmJ, wmxJ).', '%6.2f'));

%% 5. FIGURES ------------------------------------------------------------------------------
figure(40); clf; set(gcf, 'Color', 'w', 'Name', 'Minimum cycle time');
plot(mpList, tcMin, 'o-', 'LineWidth', 1.5); grid on; hold on;
yline(tc0, 'k--', 'planned cycle');
xlabel('payload [kg]'); ylabel('minimum cycle time [s]');
title('Minimum cycle time with the chosen motors and transmissions');

figure(41); clf; set(gcf, 'Color', 'w', 'Name', 'Motor utilisation vs cycle time');
for j = 1:4
    subplot(2,2,j); hold on; grid on; set(gca, 'XScale', 'log');
    for ip = [1 nP]
        plot(sList*tc0, squeeze(U(j,ip,:)), 'LineWidth', 1.3);
    end
    yline(1, 'r--', 'limit');  xline(tc0, 'k:', 'planned');
    xlabel('cycle time [s]'); ylabel('utilisation u'); ylim([0 2]);
    title(sprintf('%s - %s', jName{j}, motName{k(j)}));
    if j == 1, legend(sprintf('%.1f kg', mpList(1)), sprintf('%.1f kg', mpList(nP)), 'Location', 'best'); end
end
sgtitle('u = max(C_{pk}/C_{max}, C_{rms}/C_n, \omega/\omega_{max}) : OK while u < 1');

figure(42); clf; set(gcf, 'Color', 'w', 'Name', 'Utilisation vs payload');
u0 = zeros(4, nP);
for ip = 1:nP, u0(:,ip) = util(mpList(ip), 1, G, D, Pg, Pa, Qp, Qpp, T, tau, CnJ, CmxJ, JmJ, wmxJ); end
plot(mpList, u0, 'o-', 'LineWidth', 1.3); grid on; yline(1, 'r--', 'limit');
xlabel('payload [kg]'); ylabel('utilisation u'); legend(jName, 'Location', 'best');
title(sprintf('Motor utilisation at the planned cycle time (%.1f s)', tc0));

%% LOCAL FUNCTION ---------------------------------------------------------------------------
function [u, upk, urms, uw] = util(mp, s, G, D, Pg, Pa, Qp, Qpp, T, tau, Cn, Cmax, Jm, wmax)
% utilisation of the 4 motors for payload mp and cycle time x s
Cr   = G + D/s^2 + mp*(Pg + Pa/s^2);                    % load torque Cr* (exact scaling)
Cm   = tau.*Cr + Jm.*(Qpp/s^2)./tau;                    % course: Cm = tau*Cr* + Jm*wd_r/tau
upk  = max(abs(Cm), [], 2)./Cmax;
urms = sqrt(trapz(T, Cm.^2, 2)/(T(end) - T(1)))./Cn;    % RMS unchanged by the time scaling of the axis
uw   = (max(abs(Qp), [], 2)/s./tau)./wmax;
u    = max([upk, urms, uw], [], 2);
end