% sweep_full.m  —  Comprehensive timing + accuracy study
%
% δ ∈ {4,8,16,32}, N ∈ {4,8,16,32}, α ∈ 0:0.1:1.0
% Force model: J6dragSRPmoon   |   Span: 7 days
% All 30 orbits.
%
% Purpose: efficiency comparison — scatter of (fevals, r_err_final) for
% all MPCM δ×N×α combinations vs ode45 at 4 tolerance settings.
%
% Output: results/full/full_<name>_J6dragSRPmoon.mat  (30 files)
%   full_results   : struct array (176×1): .delta .N .alpha .gbl
%   ode45_results  : struct array (4×1):   .tol_rel .tol_abs .wall_time .fevals .r_err_final
%   delta_vec, N_vec, alpha_vec, orbit, force_model

clc;
load('benchmark_suite.mat', 'suite');
[~, mu, Re, ~, ~, ~, ~, ~] = getgravc(84);

C.mu      = mu;
C.Re      = Re;
C.J       = [0, 1.08262617385222e-3, -2.53241051856772e-6, -1.61989759991697e-6, ...
                -2.27753590730836e-7, 5.40666576283813e-7];
C.CD_drag = 2.2;
C.A_drag  = 4e-6;       % km²
C.m_drag  = 260;        % kg
C.jdepoch = 2461042.0;  % Jan 1, 2026 12:00 TT
C.rhoSRP  = 4.56e-6;    % N/m²
C.CR      = 1.8;
C.ASRP    = 4e-6;       % km²
C.muM     = 4902.8;     % km³/s²
C.debug   = false;

delta_vec  = [4, 8, 16, 32];
N_vec      = [4, 8, 16, 32];
alpha_vec  = 0 : 0.1 : 1.0;     % 11 values
force_model = 'J6dragSRPmoon';
tol_list   = [1e-6; 1e-9; 1e-13];
span_7d    = 7 * 86400;          % 604800 s

n_delta = numel(delta_vec);
n_N     = numel(N_vec);
n_alpha = numel(alpha_vec);
n_runs  = n_delta * n_N * n_alpha;  % 176 per orbit

if ~exist('results/full', 'dir'), mkdir('results/full'); end

f_ode_srp = @(t,x) orbit_eq_J6_drag_SRP_moon(t, x, C.mu, C.CD_drag, C.A_drag, C.m_drag, ...
                       C.Re, C.J, C.jdepoch, C.rhoSRP, C.CR, C.ASRP, C.muM);

for orb_idx = 1:numel(suite)
    orb = suite(orb_idx);
    x0c = [orb.r0, orb.v0]';

    %% Load 7-day slice of reference trajectory
    ref = load(fullfile('ref_trajectories', sprintf('ref_%s.mat', orb.name)));
    mask     = ref.t_ref <= span_7d;
    t_ref_7d = ref.t_ref(mask);
    X_ref_7d = ref.X_ref(mask, :);

    %% ode45 baseline (3 matched tolerances + 1 mixed)
    n_tol45       = numel(tol_list) + 1;
    ode45_results = repmat(struct('tol_rel',[], 'tol_abs',[], ...
                                  'wall_time',[], 'fevals',[], 'r_err_final',[]), ...
                           n_tol45, 1);
    for ti = 1:numel(tol_list)
        tol    = tol_list(ti);
        opts45 = odeset('RelTol', tol, 'AbsTol', tol);
        tw = inf;
        for rep = 1:3
            tic; [t45, x45] = ode45(f_ode_srp, [0, span_7d], x0c, opts45); tw = min(tw, toc);
        end
        ode45_results(ti).tol_rel     = tol;
        ode45_results(ti).tol_abs     = tol;
        ode45_results(ti).wall_time   = tw;
        ode45_results(ti).fevals      = (numel(t45) - 1) * 6;
        ode45_results(ti).r_err_final = norm(x45(end,1:3) - X_ref_7d(end,1:3));
        fprintf('[%s/J6dragSRPmoon] ode45 %.0e/%.0e  r_err=%.3e km  t=%.2f s  fevals=%d\n', ...
            orb.name, tol, tol, ode45_results(ti).r_err_final, tw, ode45_results(ti).fevals);
    end
    % Mixed 1e-6/1e-9 (MPCM-matched tolerances)
    opts_mix = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    tw = inf;
    for rep = 1:3
        tic; [t45, x45] = ode45(f_ode_srp, [0, span_7d], x0c, opts_mix); tw = min(tw, toc);
    end
    ode45_results(end).tol_rel     = 1e-6;
    ode45_results(end).tol_abs     = 1e-9;
    ode45_results(end).wall_time   = tw;
    ode45_results(end).fevals      = (numel(t45) - 1) * 6;
    ode45_results(end).r_err_final = norm(x45(end,1:3) - X_ref_7d(end,1:3));
    fprintf('[%s/J6dragSRPmoon] ode45 1e-6/1e-9  r_err=%.3e km  t=%.2f s  fevals=%d\n', ...
        orb.name, ode45_results(end).r_err_final, tw, ode45_results(end).fevals);

    %% MPCM sweep
    run_idx      = 0;
    full_results = repmat(struct('delta',[], 'N',[], 'alpha',[], 'gbl',[]), n_runs, 1);

    for di = 1:n_delta
        for ni = 1:n_N
            for ai = 1:n_alpha
                run_idx = run_idx + 1;
                fprintf('[%s/J6dragSRPmoon] d=%2d N=%2d a=%.1f  (%3d/%d) ... ', ...
                    orb.name, delta_vec(di), N_vec(ni), alpha_vec(ai), run_idx, n_runs);
                [gbl, ~, ~] = run_sweep_case(orb, force_model, alpha_vec(ai), ...
                                              delta_vec(di), N_vec(ni), ...
                                              t_ref_7d, X_ref_7d, C);
                full_results(run_idx).delta = delta_vec(di);
                full_results(run_idx).N     = N_vec(ni);
                full_results(run_idx).alpha = alpha_vec(ai);
                full_results(run_idx).gbl   = gbl;
                if gbl.error_code == 0
                    fprintf('r_err=%.3e km  t=%.2f s  fevals=%d\n', ...
                        gbl.r_err_final, gbl.wall_time, gbl.total_fevals);
                else
                    fprintf('FAIL ec=%d  t=%.2f s\n', gbl.error_code, gbl.wall_time);
                end
            end
        end
    end

    orbit = orb;
    fname = fullfile('results', 'full', sprintf('full_%s_%s.mat', orb.name, force_model));
    save(fname, 'full_results', 'ode45_results', ...
         'delta_vec', 'N_vec', 'alpha_vec', 'orbit', 'force_model');
    fprintf('  -> Saved %s\n\n', fname);
end

fprintf('sweep_full done.\n');
