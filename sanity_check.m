% sanity_check.m
%
% Runs two cases (LEO-1/J6, alpha=0 and alpha=0.5) and verifies
% run_sweep_case + odeMPGI before committing to the full sweep.
%
% All checks print PASS or FAIL. Fix any FAILs before running sweep_alpha.
clc;
load('benchmark_suite.mat', 'suite');

[~, mu1, Re, ~, ~, ~, ~, ~] = getgravc(84);
C.Re     = Re;
C.mu     = mu1;
C.J      = [0, 1.08262617385222e-3, -2.53241051856772e-6, -1.61989759991697e-6, ...
               -2.27753590730836e-7, 5.40666576283813e-7];
C.CD_drag = 2.2;
C.A_drag  = 4e-6;
C.m_drag  = 260;

delta      = 16;
N_nodes    = 16;
Sec        = suite(1).T_orbital / delta;
debug_plots = true;   % set false to skip orbital-element debug figures

ref = load('ref_trajectories/ref_LEO-1_J6.mat');

C.debug   = debug_plots;
C.debug_E = true;

fprintf('Running alpha=0  (Chebyshev branch) ...\n');
[g0, s0, n0] = run_sweep_case(suite(1), 'J6', 0,   delta, N_nodes, ref.t_ref, ref.X_ref, C);
C.debug = false;      % only plot once
fprintf('Running alpha=0.5 (Gegenbauer branch) ...\n');
[g5, s5, n5] = run_sweep_case(suite(1), 'J6', 0.5, delta, N_nodes, ref.t_ref, ref.X_ref, C);
fprintf('Running alpha=0  (Keplerian, no perturbations) ...\n');
C_kep = C; C_kep.debug = false;
[gK, sK, nK] = run_sweep_case(suite(1), 'Kep', 0, delta, N_nodes, ref.t_ref, ref.X_ref, C_kep);

%% Checks
fprintf('\n--- Sanity checks ---\n');

ok1 = (g0.error_code == 0) && (g5.error_code == 0);
fprintf('  %s  error_code==0 for both cases  (g0=%d, g5=%d)\n', ...
    tf(ok1), g0.error_code, g5.error_code);

ok2 = isfinite(g0.r_err_final) && g0.r_err_final > 0 && ...
      isfinite(g5.r_err_final) && g5.r_err_final > 0;
fprintf('  %s  r_err_final finite and positive  (g0=%.3e km, g5=%.3e km)\n', ...
    tf(ok2), g0.r_err_final, g5.r_err_final);

ok3 = (sum(s0.iters) == g0.total_iters) && (sum(s5.iters) == g5.total_iters);
fprintf('  %s  sum(seg.iters) == total_iters  (g0: %d==%d, g5: %d==%d)\n', ...
    tf(ok3), sum(s0.iters), g0.total_iters, sum(s5.iters), g5.total_iters);

ok4 = numel(s0.r_err) == numel(s5.r_err);
fprintf('  %s  same segment count for both cases  (%d segments)\n', ...
    tf(ok4), numel(s0.r_err));

ok5 = numel(n0.r_err) == numel(ref.t_ref);
fprintf('  %s  per-node count == length(t_ref)  (%d==%d)\n', ...
    tf(ok5), numel(n0.r_err), numel(ref.t_ref));

ok6 = all(isfinite(s0.e_osc)) && all(isfinite(s0.E_osc));
fprintf('  %s  e_osc and E_osc all finite  (mean e_osc=%.4f, expect ~0.001)\n', ...
    tf(ok6), mean(s0.e_osc));

ok7 = ~any(isnan(s0.r_start)) && ~any(isnan(s5.r_start));
fprintf('  %s  no NaN in r_start\n', tf(ok7));

% ok8: segment count matches expected value from t_ref span (±1 for boundary)
M_expected = floor(ref.t_ref(end) / Sec);
M0 = numel(s0.r_err);
M5 = numel(s5.r_err);
ok8 = (M0 == M5) && abs(M0 - M_expected) <= 1;
fprintf('  %s  seg count matches t_ref span  (expected~%d, s0=%d, s5=%d)\n', ...
    tf(ok8), M_expected, M0, M5);

% ok9: all seg fields finite (no NaN or Inf anywhere)
seg_fields = {'r_err','v_err','iters','fevals','t_start','r_start','v_start','e_osc','f_osc','E_osc'};
ok9 = all(cellfun(@(f) all(isfinite(s0.(f))), seg_fields)) && ...
      all(cellfun(@(f) all(isfinite(s5.(f))), seg_fields));
fprintf('  %s  all seg fields finite (no NaN/Inf)\n', tf(ok9));

% ok10: seg.t_start strictly increasing
ok10 = all(diff(s0.t_start) > 0) && all(diff(s5.t_start) > 0);
fprintf('  %s  seg.t_start strictly increasing\n', tf(ok10));

% ok11: seg.t_start spacing == Sec within 0.1%
tol = 1e-3 * Sec;
dev0 = max(abs(diff(s0.t_start) - Sec));
dev5 = max(abs(diff(s5.t_start) - Sec));
ok11 = (dev0 < tol) && (dev5 < tol);
fprintf('  %s  seg.t_start spacing ~= Sec=%.3fs  (max dev: s0=%.2e s, s5=%.2e s)\n', ...
    tf(ok11), Sec, dev0, dev5);

% ok12: raw E_osc — every non-positive diff must be an atan2 wraparound
%       (diff < -pi, i.e. near -2pi). Any decrease that is NOT a wraparound
%       means E_osc is genuinely non-monotonic and unwrap cannot fix it.
d0raw = diff(s0.E_osc);
d5raw = diff(s5.E_osc);
ok12 = all((d0raw > 0) | (d0raw < -pi)) && ...
       all((d5raw > 0) | (d5raw < -pi));
bad0 = find(d0raw <= 0 & d0raw >= -pi);
bad5 = find(d5raw <= 0 & d5raw >= -pi);
fprintf('  %s  raw E_osc jumps only at atan2 boundary (+/-pi)  (bad segs s0: [%s], s5: [%s])\n', ...
    tf(ok12), num2str(bad0(:)'), num2str(bad5(:)'));

% ok13 (was ok12): E_osc monotonically increasing after unwrap
E0u = unwrap(s0.E_osc);
E5u = unwrap(s5.E_osc);
ok13_E = all(diff(E0u) > 0) && all(diff(E5u) > 0);
fprintf('  %s  E_osc monotonically increasing after unwrap\n', tf(ok13_E));

% ok13: seg.fevals == seg.iters * N_nodes element-wise
ok13 = all(s0.fevals == s0.iters * N_nodes) && ...
       all(s5.fevals == s5.iters * N_nodes);
fprintf('  %s  seg.fevals == seg.iters * N_nodes (element-wise)\n', tf(ok13));

% ok14: gbl.total_fevals == sum(seg.iters) * N_nodes
ok14 = (g0.total_fevals == sum(s0.iters) * N_nodes) && ...
       (g5.total_fevals == sum(s5.iters) * N_nodes);
fprintf('  %s  total_fevals == sum(seg.iters)*N  (g0: %d==%d, g5: %d==%d)\n', ...
    tf(ok14), g0.total_fevals, sum(s0.iters)*N_nodes, ...
              g5.total_fevals, sum(s5.iters)*N_nodes);

% ok15: no NaN/Inf in nod fields
ok15 = all(isfinite(n0.r_err)) && all(isfinite(n0.v_err)) && ...
       all(isfinite(n5.r_err)) && all(isfinite(n5.v_err));
fprintf('  %s  no NaN/Inf in nod.r_err / nod.v_err\n', tf(ok15));

% ok16: nod.t strictly increasing (no repeated or out-of-order time points)
ok16 = all(diff(n0.t) > 0) && all(diff(n5.t) > 0);
fprintf('  %s  nod.t strictly increasing (no repeated points)\n', tf(ok16));

% ok17: nod.t matches t_ref exactly
ok17 = isequal(n0.t, ref.t_ref(:)) && isequal(n5.t, ref.t_ref(:));
fprintf('  %s  nod.t == t_ref\n', tf(ok17));

% ok18: nod.seg_idx all in [1, M] (M = actual segment count)
ok18 = all(n0.seg_idx >= 1 & n0.seg_idx <= M0) && ...
       all(n5.seg_idx >= 1 & n5.seg_idx <= M5);
fprintf('  %s  nod.seg_idx all in [1,M]  (s0 max=%d vs M=%d, s5 max=%d vs M=%d)\n', ...
    tf(ok18), max(n0.seg_idx), M0, max(n5.seg_idx), M5);

% ok19: every segment 1:M has at least one node
segs0_present = numel(unique(n0.seg_idx(n0.seg_idx >= 1 & n0.seg_idx <= M0)));
segs5_present = numel(unique(n5.seg_idx(n5.seg_idx >= 1 & n5.seg_idx <= M5)));
ok19 = (segs0_present == M0) && (segs5_present == M5);
fprintf('  %s  all M segments have >= 1 node  (s0=%d/%d, s5=%d/%d covered)\n', ...
    tf(ok19), segs0_present, M0, segs5_present, M5);

%% Diagnostic tables — E1/E2/E3 comparison, first 2 orbits only
fprintf('\n--- E comparison: J6 alpha=0 (segs 1..%d) ---\n', 2*delta);
fprintf('%4s  %10s  %10s  %10s  %10s  %10s\n','seg','t_start','E1','E2','E3','r_err');
for k = 1:min(2*delta, numel(s0.t_start))
    fprintf('%4d  %10.3f  %10.6f  %10.6f  %10.6f  %10.3e\n', ...
        k, s0.t_start(k), s0.E1_osc(k), s0.E2_osc(k), s0.E3_osc(k), s0.r_err(k));
end

fprintf('\n--- E comparison: Keplerian alpha=0 (segs 1..%d) ---\n', 2*delta);
fprintf('%4s  %10s  %10s  %10s  %10s\n','seg','t_start','E1','E2','E3');
for k = 1:min(2*delta, numel(sK.t_start))
    fprintf('%4d  %10.3f  %10.6f  %10.6f  %10.6f\n', ...
        k, sK.t_start(k), sK.E1_osc(k), sK.E2_osc(k), sK.E3_osc(k));
end

%% E comparison figure
t_orb = (0:2*delta-1)' / delta;
figure('Name','E comparison: Keplerian vs J6','Position',[50 50 1200 800]);
labels = {'E1 (orbit-eq)','E2 (atan2)','E3 (Kepler eq)'};
flds   = {'E1_osc','E2_osc','E3_osc'};
runs   = {sK, s0};
ttls   = {'Keplerian','J6'};
for col = 1:2
    sr = runs{col};
    for row = 1:3
        subplot(3,2,(row-1)*2+col);
        n_show = min(2*delta, numel(sr.(flds{row})));
        plot(t_orb(1:n_show), rad2deg(sr.(flds{row})(1:n_show)), '.-');
        xlabel('Orbit fraction'); ylabel('E (deg)');
        title(sprintf('%s — %s', labels{row}, ttls{col}));
        grid on;
    end
end
sgtitle('E formula comparison (first 2 orbits)');

fprintf('\n--- kep_elements spot-check (alpha=0, segs 1,4,8,12,16) ---\n');
fprintf('%4s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n', ...
    'seg', 'a', 'e', 'i_deg', 'w_deg', 'f_deg', 'E_deg', 'u_deg');
for k = [1, 4, 8, 12, 16]
    [ak,ek,ik,~,wk,fk] = kep_elements(s0.rv_start(k,:), s0.vv_start(k,:), C.mu);
    Ek = atan2(sqrt(max(0,1-ek^2))*sin(fk), ek+cos(fk));
    uk = mod(wk + fk, 2*pi);
    fprintf('%4d  %8.2f  %8.5f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n', ...
        k, ak, ek, rad2deg(ik), rad2deg(wk), rad2deg(fk), rad2deg(Ek), rad2deg(uk));
end

fprintf('\n--- nod diagnostics (alpha=0, first and last 5 points) ---\n');
fprintf('%6s  %12s  %8s  %10s\n', 'idx', 't', 'seg_idx', 'r_err');
show = unique([1:5, max(1,numel(n0.t)-4):numel(n0.t)]);
for k = show
    fprintf('%6d  %12.4f  %8d  %10.3e\n', k, n0.t(k), n0.seg_idx(k), n0.r_err(k));
end

%% Summary
checks = [ok1,ok2,ok3,ok4,ok5,ok6,ok7,ok8,ok9,ok10,ok11,ok12,ok13_E,ok13,ok14,ok15,ok16,ok17,ok18,ok19];
fprintf('\n%d / %d checks passed.\n', sum(checks), numel(checks));
if all(checks)
    fprintf('All checks passed — safe to run sweep_alpha.\n');
else
    fprintf('Fix failing checks before running sweep_alpha.\n');
end

function s = tf(ok)
    if ok; s = 'PASS'; else; s = 'FAIL'; end
end
