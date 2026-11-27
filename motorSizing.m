% motorSizing.m - motor + reducer choice with the "dynamic sizing" of the course
% (Giberti, Motor Sizing Notes). Run AFTER projectA_main (uses T, Q, Qp, Qpp, atip, L_dim, alpha, Masses).
% Course notation: omega_r = tau*omega_m (tau < 1 = speed reduction)
%   Cm/tau - Cr = (Jm/tau^2 + Jr)*wd_r      ->   Cm = tau*Cr* + Jm*wd_r/tau          (slides 9, 21)
%   load torque     Cr* = Cr + Jr*wd_r      = joint effort from the dynamic model (+ payload)
%   load factor     beta  = 2*[wd_rq*C*_rq + mean(wd_r*C*_r)]                     (slide 24)
%   accel. factor   alpha = Cn^2/Jm   -> the motor can be used only if alpha > beta (slides 24, 26)
%   tau_min,max = sqrt(Jm)*(sqrt(alpha-beta+4*wd_rq*C*_rq) -/+ sqrt(alpha-beta))/(2*C*_rq)  (slide 26)
%   tau_opt = sqrt(Jm*wd_rq/C*_rq)                                                  (slide 25)
%   tau_p   = w_r,max/w_m,max   kinematic limit, tau > tau_p                        (slides 16, 27)
%   checks: Cmax > max|tau*Cr* + Jm*wd_r/tau|,  Crms < Cn,  w_max < w_m,max         (slides 10, 19, 28)
%   reducer: maximum input velocity, maximum output torque                          (slide 10)
% q2 (prismatic): same formulas with mass/force/linear speed (slide 32); tau = lead/(2*pi) [m/rad].
if ~exist('Masses', 'var') || ~exist('atip', 'var'), error('Run projectA_main first.'); end

%% 1. DATA ------------------------------------------------------------------
% allowed motors (project sheet: HG-MR / HG-KR, <= 200 W). Order = preference at equal power
motName = {'HG-KR053', 'HG-MR053', 'HG-KR13', 'HG-MR13', 'HG-KR23', 'HG-MR23'};
Cn      = [0.16   0.16   0.32   0.32   0.64   0.64  ];      % rated torque   [N m]
CmMax   = [0.56   0.48   1.1    0.95   2.2    1.9   ];      % maximum torque [N m]
Jm      = [0.0450 0.0162 0.0777 0.0300 0.221  0.0865]*1e-4; % rotor inertia  [kg m^2] (no brake)
wmMax   = 6000*2*pi/60*ones(1, 6);                           % maximum speed  [rad/s]
nMot    = numel(motName);

payloadSize  = 1.0;             % sizing with the worst payload of the sheet [kg]
payloadCheck = 0.2;             % then check with the lightest one [kg]
margin    = 1.2;                % [ASSUMPTION] safety margin: tau = (largest admissible tau)/margin
iMaxPract = 60;                 % [ASSUMPTION] largest reduction 1/tau you can build (rotary joints)
leadLim   = [0.002 0.020];      % [ASSUMPTION] screw leads you can buy/build for q2 [m]
jName = {'q1 (base)', 'q2 (prismatic)', 'q3 (elbow)', 'q4 (wrist)'};
tc = T(end) - T(1);                                         % cycle time

%% 2. LOAD SIDE: Cr* along the task (dynamic model + payload at the TCP) ------
F0  = InverseDynamics(Q(:,1), zeros(4,1), zeros(4,1), L_dim, alpha, Masses);
sgn = sign(F0(2));  if sgn == 0, sgn = 1; end               % same sign convention as InverseDynamics
CrS = loadEfforts(Q, Qp, Qpp, atip, L_dim, alpha, Masses, payloadSize,  sgn);   % 4 x M
CrC = loadEfforts(Q, Qp, Qpp, atip, L_dim, alpha, Masses, payloadCheck, sgn);

%% 3. SIZING JOINT BY JOINT -----------------------------------------------------
S = struct([]);
for j = 1:4
    Cr = CrS(j,:);  wd = Qpp(j,:);  w = Qp(j,:);
    Crq   = sqrt(trapz(T, Cr.^2)/tc);                       % C*_r,q
    wdq   = sqrt(trapz(T, wd.^2)/tc);                       % wd_r,q
    CwM   = trapz(T, wd.*Cr)/tc;                            % (wd_r C*_r) mean
    beta  = 2*(wdq*Crq + CwM);                              % load factor [W/s]
    wrMax = max(abs(w));
    if j == 2, tPr = leadLim/(2*pi); else, tPr = [1/iMaxPract, 1]; end   % practical tau range

    r = struct('Crq', Crq, 'wdq', wdq, 'beta', beta, 'wrMax', wrMax, 'alpha', Cn.^2./Jm, ...
               'tmin', nan(1,nMot), 'tmax', nan(1,nMot), 'topt', nan(1,nMot), 'tp', nan(1,nMot), ...
               'tau', nan(1,nMot), 'Cpk', nan(1,nMot), 'Crms', nan(1,nMot), 'wmPk', nan(1,nMot), ...
               'ok', false(1,nMot), 'tPr', tPr, 'pick', 0);
    for k = 1:nMot
        r.tp(k)   = wrMax/wmMax(k);
        r.topt(k) = sqrt(Jm(k)*wdq/Crq);
        d = r.alpha(k) - beta;
        if d <= 0, continue; end                            % alpha <= beta: motor too weak
        r.tmin(k) = sqrt(Jm(k))*(sqrt(d + 4*wdq*Crq) - sqrt(d))/(2*Crq);
        r.tmax(k) = sqrt(Jm(k))*(sqrt(d + 4*wdq*Crq) + sqrt(d))/(2*Crq);
        lo = max([r.tmin(k), r.tp(k), tPr(1)]);  hi = min(r.tmax(k), tPr(2));
        if lo >= hi, continue; end
        tg  = logspace(log10(lo), log10(hi), 400);         % admissible tau: also the peak check
        pk  = arrayfun(@(t) max(abs(t*Cr + Jm(k)*wd/t)), tg);
        okT = pk < CmMax(k);
        if ~any(okT), continue; end
        tHi = max(tg(okT));  tLo = min(tg(okT));
        r.tau(k) = max(tLo, tHi/margin);                    % largest tau (least reduction) with margin
        Cm = r.tau(k)*Cr + Jm(k)*wd/r.tau(k);               % motor torque along the task
        r.Cpk(k)  = max(abs(Cm));
        r.Crms(k) = sqrt(trapz(T, Cm.^2)/tc);
        r.wmPk(k) = wrMax/r.tau(k);
        r.ok(k)   = r.Cpk(k) < CmMax(k) && r.Crms(k) < Cn(k) && r.wmPk(k) < wmMax(k);
    end
    p = find(r.ok, 1);                                      % smallest suitable motor
    if ~isempty(p), r.pick = p; end
    if j == 1, S = r; else, S(j) = r; end
end

%% 4. PRINT ------------------------------------------------------------------------
fprintf('\n=== LOAD SIDE (payload %.1f kg, cycle %.1f s) ===\n', payloadSize, tc);
for j = 1:4
    fprintf('%-15s C*rq = %8.3f   wd_rq = %8.4f   beta = %9.4f W/s   w_r,max = %.4f\n', ...
            jName{j}, S(j).Crq, S(j).wdq, S(j).beta, S(j).wrMax);
end
fprintf('(q2: C* in N, wd in m/s^2, w in m/s; tau in m/rad = lead/2pi)\n');

for j = 1:4
    fprintf('\n--- %s ---\n%-9s %10s %6s %10s %10s %10s %10s %10s %8s %8s %8s  %s\n', jName{j}, ...
            'motor', 'alpha', '>beta', 'tau_min', 'tau_max', 'tau_opt', 'tau_p', 'tau', 'ratio', ...
            'Cpk', 'Crms', 'checks');
    for k = 1:nMot
        if j == 2, rt = sprintf('%6.1fmm', 1e3*2*pi*S(j).tau(k)); else, rt = sprintf('%8.1f', 1/S(j).tau(k)); end
        fprintf('%-9s %10.3g %6s %10.3g %10.3g %10.3g %10.3g %10.3g %8s %8.3f %8.3f  %s\n', motName{k}, ...
                S(j).alpha(k), yn(S(j).alpha(k) > S(j).beta), S(j).tmin(k), S(j).tmax(k), S(j).topt(k), ...
                S(j).tp(k), S(j).tau(k), rt, S(j).Cpk(k), S(j).Crms(k), okTxt(S(j).ok(k)));
    end
end

fprintf('\n=== CHOICE (smallest suitable motor, tau with margin %.2f) + check with %.1f kg ===\n', margin, payloadCheck);
for j = 1:4
    k = S(j).pick;
    if k == 0, fprintf('%-15s NO SUITABLE MOTOR: raise iMaxPract / leadLim or lighten the robot\n', jName{j}); continue; end
    t  = S(j).tau(k);
    Cm = t*CrC(j,:) + Jm(k)*Qpp(j,:)/t;
    okC = max(abs(Cm)) < CmMax(k) && sqrt(trapz(T, Cm.^2)/tc) < Cn(k);
    if j == 2, red = sprintf('screw lead %.1f mm', 1e3*2*pi*t); else, red = sprintf('reduction 1/tau = %.1f', 1/t); end
    fprintf('%-15s %-9s %-24s | reducer: input <= %5.0f rpm, output <= %8.2f %s | %.1f kg: %s\n', ...
            jName{j}, motName{k}, red, S(j).wrMax/t*60/(2*pi), max(abs(CrS(j,:))), unitTxt(j), ...
            payloadCheck, okTxt(okC));
end

%% 5. FIGURES (as slides 29, 33, 34) -------------------------------------------------
figure(30); clf; set(gcf, 'Color', 'w', 'Name', 'alpha-beta diagram');
for j = 1:4
    subplot(2,2,j); hold on; grid on; set(gca, 'YScale', 'log');
    plot(1:nMot, S(j).alpha, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.35 0.6 0.9]);
    yline(S(j).beta, 'r-', '\beta', 'LineWidth', 1.5);
    if S(j).pick > 0, plot(S(j).pick, S(j).alpha(S(j).pick), 'ro', 'MarkerFaceColor', 'r'); end
    set(gca, 'XTick', 1:nMot, 'XTickLabel', motName, 'XTickLabelRotation', 30);
    ylabel('\alpha = C_n^2/J_m  [W/s]'); title(jName{j});
end
sgtitle('\alpha-\beta diagram (red dot: chosen motor)');

figure(31); clf; set(gcf, 'Color', 'w', 'Name', 'tau diagram');
for j = 1:4
    subplot(2,2,j); hold on; grid on; set(gca, 'YScale', 'log');
    h = [plot(1:nMot, S(j).tmax, 'b^'), plot(1:nMot, S(j).tmin, 'bv'), plot(1:nMot, S(j).topt, 'go'), ...
         plot(1:nMot, S(j).tp, 'k*'), plot(1:nMot, S(j).tau, 'rs', 'MarkerFaceColor', 'r')];
    yline(S(j).tPr(1), 'm--');  yline(S(j).tPr(2), 'm--');
    set(gca, 'XTick', 1:nMot, 'XTickLabel', motName, 'XTickLabelRotation', 30);
    ylabel('\tau'); title(jName{j});
    if j == 1, legend(h, {'\tau_{max}', '\tau_{min}', '\tau_{opt}', '\tau_p', '\tau chosen'}, 'Location', 'best'); end
end
sgtitle('\tau diagram (magenta: practical range of your reducer)');

figure(32); clf; set(gcf, 'Color', 'w', 'Name', 'Chosen motors: torque-speed');
for j = 1:4
    subplot(2,2,j); hold on; grid on;
    k = S(j).pick;  if k == 0, title([jName{j} ': no motor']); continue; end
    t  = S(j).tau(k);  Cm = t*CrS(j,:) + Jm(k)*Qpp(j,:)/t;  wm = Qp(j,:)/t*60/(2*pi);
    nM = wmMax(k)*60/(2*pi);
    plot(wm, Cm, 'b', 'LineWidth', 1.2);
    plot([-nM nM], [1 1]*Cn(k), 'g', [-nM nM], -[1 1]*Cn(k), 'g', 'LineWidth', 1.5);
    plot([-nM nM], [1 1]*CmMax(k), 'r', [-nM nM], -[1 1]*CmMax(k), 'r', 'LineWidth', 1.5);
    yline(S(j).Crms(k), 'b--', 'C_{rms}');  xline(nM, 'r--');  xline(-nM, 'r--');
    xlabel('\omega_m [rpm]'); ylabel('C_m [N m]');
    title(sprintf('%s: %s, %s', jName{j}, motName{k}, ratioTxt(j, t)));
end
sgtitle('Chosen motors: C_n (green), C_{max} (red), task (blue)');

%% LOCAL FUNCTIONS ---------------------------------------------------------------------
function Cr = loadEfforts(Q, Qp, Qpp, atip, L_dim, alpha, Masses, mp, sgn)
% Cr* = dynamic model (links) + payload mp at the TCP (J' * force on the payload)
M = size(Q, 2);  Cr = zeros(4, M);
for k = 1:M
    J  = Robot4DOF_jac(Q(:,k), L_dim, alpha);
    Fr = InverseDynamics(Q(:,k), Qp(:,k), Qpp(:,k), L_dim, alpha, Masses);
    Fp = mp*(atip(:,k) + [0; 0; 9.81]);
    Cr(:,k) = Fr(:) + sgn*J(1:3,:).'*Fp;
end
end
function s = yn(b),    if b, s = 'yes'; else, s = 'NO'; end, end
function s = okTxt(b), if b, s = 'OK'; else, s = 'not OK'; end, end
function s = unitTxt(j), if j == 2, s = 'N'; else, s = 'N m'; end, end
function s = ratioTxt(j, t)
if j == 2, s = sprintf('lead %.1f mm', 1e3*2*pi*t); else, s = sprintf('1/\\tau = %.1f', 1/t); end
end