% test_sweep_alpha.m
%
% Smoke test: builds the HEO-E07 reference trajectory exactly as
% init_benchmark_suite does, then runs sweep_alpha on HEO-E07 at
% alpha in {0, 0.5, 1.0} (delta=16, N=16).  Takes ~5-15 min.
%
% Writes to the same paths as the full pipeline so subsequent plot scripts
% work on the single-orbit result without modification.
%
% PASS criteria:
%   ref  — |r_final| in plausible LEO range, no NaN/Inf
%   runs — all 3 sweep cases return error_code=0 and finite r_err < 100 km

clc;
[~, mu, Re, ~, ~, ~, ~, ~] = getgravc(84);
J = [0, 1.08262617385222e-3, -2.53241051856772e-6, -1.61989759991697e-6, ...
        -2.27753590730836e-7, 5.40666576283813e-7];

C.mu      = mu;
C.Re      = Re;
C.J       = J;
C.CD_drag = 2.2;
C.A_drag  = 4e-6;       % km²
C.m_drag  = 260;        % kg
C.jdepoch = 2461042.0;  % Jan 1, 2026 12:00 TT
C.rhoSRP  = 4.56e-6;    % N/m²
C.CR      = 1.8;
C.ASRP    = 4e-6;       % km²
C.muM     = 4902.8;     % km³/s²

%% Build HEO-E07 orbit struct (matches init_benchmark_suite exactly)
% perigee=500 km, apogee=35786 km, i=28.5 deg
hp  = 500;   ha = 35786;   inc = 28.5;
rp  = Re + hp;
ra  = Re + ha;
a   = (rp + ra) / 2;
e   = (ra - rp) / (ra + rp);
[r0, v0] = posnvelos(a, e, inc*pi/180, 0, 0, 0, mu);
orb.name      = 'HEO-E07';
orb.r0        = r0;
orb.v0        = v0;
orb.T_orbital = 2*pi * sqrt(a^3 / mu);
orb.e         = e;
orb.h_perigee = hp;
orb.inc       = inc;
orb.sweep     = 'ecc';

fprintf('HEO-E07: a=%.1f km  e=%.4f  i=%.1f deg  T=%.1f s\n', ...
    a, e, inc, orb.T_orbital);

%% Quick force evaluation at r0 — checks the function runs without error
f_ode = @(t,x) orbit_eq_J6_drag_SRP_moon(t, x, mu, C.CD_drag, C.A_drag, C.m_drag, ...
                   Re, J, C.jdepoch, C.rhoSRP, C.CR, C.ASRP, C.muM);
dr0 = f_ode(0, [r0, v0]');
a_norm = norm(dr0(4:6));
fprintf('Force eval at t=0: |a| = %.4e km/s²  (two-body ~%.4e)\n', ...
    a_norm, mu/a^2);

%% Generate reference trajectory (same settings as init_benchmark_suite)
t_span = 30 * 86400;
dt     = orb.T_orbital / 2048;
t_eval = (0 : dt : t_span)';
if t_eval(end) < t_span, t_eval(end+1) = t_span; end

opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
fprintf('\nGenerating HEO-E07 reference trajectory (30 days, ode78) ...\n');
tic;
[t_ref, X_ref] = ode78(f_ode, t_eval, [r0, v0], opts);
elapsed = toc;
fprintf('Done in %.1f s  (%d steps)\n', elapsed, numel(t_ref));

r_final = norm(X_ref(end, 1:3));
ok_ref  = isfinite(r_final) && r_final > rp*0.9 && r_final < ra*1.1;
fprintf('|r_final| = %.3f km  [expected ~%.0f km]  %s\n', ...
    r_final, a, pass_fail(ok_ref));

%% Plot reference trajectory with Earth
figure('Name', 'HEO-E07 reference trajectory', 'Position', [50 50 700 700]);
ax = axes;
hold(ax, 'on'); axis(ax, 'equal'); grid(ax, 'on');

% Earth sphere
[sx, sy, sz] = sphere(64);
surf(ax, sx*Re, sy*Re, sz*Re, ...
    'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

% Orbit
plot3(ax, X_ref(:,1), X_ref(:,2), X_ref(:,3), 'w-', 'LineWidth', 0.8);

xlabel(ax, 'X (km)'); ylabel(ax, 'Y (km)'); zlabel(ax, 'Z (km)');
title(ax, sprintf('HEO-E07 — 30-day reference trajectory (e=%.3f)', e));
set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
set(gcf, 'Color', 'k');
view(ax, 30, 20);

%% Osculating elements from reference trajectory
r_all    = X_ref(:,1:3);
v_all    = X_ref(:,4:6);
r_norm_all = sqrt(sum(r_all.^2, 2));
v_norm_all = sqrt(sum(v_all.^2, 2));
t_days   = t_ref / 86400;

% Semi-major axis (vis-viva)
a_osc = 1 ./ (2./r_norm_all - v_norm_all.^2 ./ mu);

% Eccentricity vector and magnitude
h_all    = cross(r_all, v_all, 2);
e_vec_all = cross(v_all, h_all, 2) ./ mu - r_all ./ r_norm_all;
e_osc    = sqrt(sum(e_vec_all.^2, 2));

% Eccentric anomaly via Kepler's equation (E3 — same as run_sweep_case)
n_mean_ref = 2*pi / orb.T_orbital;
M_all      = mod(n_mean_ref .* t_ref, 2*pi);
E_iter     = M_all + e_osc .* sin(M_all);
for kk = 1:10
    dE     = (M_all - E_iter + e_osc .* sin(E_iter)) ./ (1 - e_osc .* cos(E_iter));
    E_iter = E_iter + dE;
    if max(abs(dE)) < 1e-14, break; end
end
E_osc = mod(E_iter + pi, 2*pi) - pi;

n_neg = sum(a_osc <= 0);
if n_neg > 0
    fprintf('WARNING: %d / %d points have a_osc <= 0.\n', n_neg, numel(a_osc));
else
    fprintf('All a_osc > 0.\n');
end

% Plot
figure('Name', 'Osculating elements — HEO-E07', 'Position', [800 50 800 600]);

subplot(3, 1, 1);
plot(t_days, a_osc, 'b-', 'LineWidth', 0.5);
ylabel('a_{osc} (km)'); xlabel('Time (days)');
title('Osculating semi-major axis'); grid on;

subplot(3, 1, 2);
plot(t_days, e_osc, 'r-', 'LineWidth', 0.5);
ylabel('e_{osc}'); xlabel('Time (days)');
title('Osculating eccentricity'); grid on;

subplot(3, 1, 3);
plot(t_days, unwrap(E_osc), 'm-', 'LineWidth', 0.5);
ylabel('E_{osc} (rad, unwrapped)'); xlabel('Time (days)');
title('Osculating eccentric anomaly'); grid on;

sgtitle('HEO-E07 osculating elements — 30-day reference');

%% Save ref file (same path as init_benchmark_suite)
if ~exist('ref_trajectories', 'dir'), mkdir('ref_trajectories'); end
orbit       = orb;
force_model = 'J6dragSRPmoon';
save(fullfile('ref_trajectories', 'ref_HEO-E07.mat'), ...
    't_ref', 'X_ref', 'orbit', 'force_model');
fprintf('Saved ref_trajectories/ref_HEO-E07.mat\n');

%% Run sweep_alpha for 3 representative alpha values
alpha_test = [0; 0.5; 1.0];
delta      = 16;
N_nodes    = 16;
fm         = 'J6dragSRPmoon';
n_alpha    = numel(alpha_test);

if ~exist('results/sweep', 'dir'), mkdir('results/sweep'); end

global_res = repmat(struct('alpha',[], 'r_err_final',[], 'v_err_final',[], ...
                           'total_iters',[], 'total_fevals',[], ...
                           'wall_time',[], 'error_code',[]), n_alpha, 1);
seg_res = cell(n_alpha, 1);
nod_res = cell(n_alpha, 1);

fprintf('\nRunning sweep_alpha (alpha in {0, 0.5, 1.0}) ...\n');
for ai = 1:n_alpha
    alpha = alpha_test(ai);
    fprintf('[HEO-E07/%s] alpha=%.1f ... ', fm, alpha);
    [global_res(ai), seg_res{ai}, nod_res{ai}] = ...
        run_sweep_case(orb, fm, alpha, delta, N_nodes, t_ref, X_ref, C);
    fprintf('r_err=%.3e km  fevals=%d  ec=%d\n', ...
        global_res(ai).r_err_final, global_res(ai).total_fevals, ...
        global_res(ai).error_code);
end

%% Save (same path/format as sweep_alpha so plot_alpha_sweep can read it)
alpha_vec = alpha_test;
orbit     = orb;
save(fullfile('results/sweep', sprintf('sweep_HEO-E07_%s.mat', fm)), ...
    'alpha_vec', 'global_res', 'seg_res', 'nod_res', ...
    'orbit', 'force_model', 'delta', 'N_nodes');
fprintf('Saved results/sweep/sweep_HEO-E07_%s.mat\n', fm);

%% PASS / FAIL summary
ok_ec   = all([global_res.error_code] == 0);
ok_err  = all(isfinite([global_res.r_err_final])) && ...
          all([global_res.r_err_final] < 100);
ok_runs = ok_ec && ok_err;

fprintf('\n--- Smoke test results ---\n');
fprintf('  ref trajectory   : %s  (|r_final|=%.1f km)\n', pass_fail(ok_ref), r_final);
fprintf('  sweep_alpha runs : %s\n', pass_fail(ok_runs));
for ai = 1:n_alpha
    fprintf('    alpha=%.1f  ec=%d  r_err=%.3e km\n', ...
        alpha_test(ai), global_res(ai).error_code, global_res(ai).r_err_final);
end

if ok_ref && ok_runs
    fprintf('\nALL PASS\n');
else
    fprintf('\nFAILURES DETECTED — check output above.\n');
end

% -------------------------------------------------------------------------
function s = pass_fail(x)
    if x, s = 'PASS'; else, s = 'FAIL'; end
end
