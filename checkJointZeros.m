% checkJointZeros.m - run AFTER projectA_main (uses L_dim, alpha, Q_home, Q_off, Q_sgn)
% Compares the Matlab FK and the Simscape model at FIXED joint configurations,
% one joint moved at a time, to find the zero offsets (Q_off) and signs (Q_sgn).
mdl = 'simulink_multibody_SCARA';
load_system(mdl);
try, set_param(mdl, 'EnablePacing', 'off'); end           % no real time for these short runs

d     = [30*pi/180; 0.05; 30*pi/180; 30*pi/180];          % step applied to each joint
cfg   = [Q_off, repmat(Q_home, 1, 5)];                    % 1: Simscape joints all at 0, 2: home
for j = 1:4, cfg(j, j+2) = cfg(j, j+2) + d(j); end        % 3..6: home with joint j moved
names = {'Simscape 0', 'home', 'q1 +30deg', 'q2 +50mm', 'q3 +30deg', 'q4 +30deg'};

nC = size(cfg, 2);  pM = zeros(3, nC);  pS = zeros(3, nC);
for c = 1:nC
    Qc  = cfg(:,c);
    M04 = ForwardKinematics(Qc, zeros(4,1), zeros(4,1), L_dim, alpha);
    pM(:,c) = M04(1:3,4);                                 % Matlab tip

    Qs     = (Qc - Q_off).*Q_sgn;                         % same conversion as the main
    Q_ts   = timeseries([Qs Qs].', [0; 0.2]);             % constant joint inputs
    Qp_ts  = timeseries(zeros(2,4), [0; 0.2]);
    Qpp_ts = timeseries(zeros(2,4), [0; 0.2]);
    o = sim(mdl, 'StopTime', '0.2');
    s = o.get('tipSim');
    if isa(s, 'timeseries'), s = s.Data; end
    s = squeeze(s);  if size(s,2) ~= 3, s = s.'; end
    pS(:,c) = s(end,:).';                                 % Simscape tip
end

fprintf('\n%-11s %-24s %-24s %s\n', 'config', 'Matlab tip [m]', 'Simscape tip [m]', 'error [mm]');
for c = 1:nC
    fprintf('%-11s [%6.3f %6.3f %6.3f]  [%6.3f %6.3f %6.3f]  %7.1f\n', ...
        names{c}, pM(:,c), pS(:,c), 1e3*norm(pM(:,c) - pS(:,c)));
end
fprintf('\nTip displacement caused by each joint (from home):\n');
for c = 3:nC
    fprintf('%-11s Matlab [%6.3f %6.3f %6.3f]   Simscape [%6.3f %6.3f %6.3f]\n', ...
        names{c}, pM(:,c) - pM(:,2), pS(:,c) - pS(:,2));
end