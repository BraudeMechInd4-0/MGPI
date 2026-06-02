% sweep_sensitivity.m  —  Stage 5.3
%
% Sensitivity sweep: δ ∈ {4,8,16,32}, N ∈ {8,16,24,32}, α ∈ {0,0.25,0.5,0.75,1.0}
% on 4 representative orbits: LEO-1, MEO-1, HEO-E07, GEO-1.
% Single force model: J6dragSRPmoon.
% Only global metrics stored (no seg/nod) — goal is (δ,N,α) selection, not per-segment analysis.
%
% Output: results/sensitivity/sens_<name>_J6dragSRPmoon.mat per rep orbit
%   sens_results : struct array — fields: delta, N, alpha, gbl
%   delta_vec, N_vec, alpha_vec5 : swept parameter vectors
%   orbit, force_model

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

delta_vec  = [4, 8, 16, 32];
N_vec      = [8, 16, 24, 32];
alpha_vec5 = [0, 0.25, 0.5, 0.75, 1.0];
force_model = 'J6dragSRPmoon';
rep_names   = {'HEO-E07', 'LEO-1', 'MEO-1', 'GEO-1'};

% Map rep_names to suite indices
rep_idx = zeros(1, numel(rep_names));
for j = 1:numel(rep_names)
    for ki = 1:numel(suite)
        if strcmp(suite(ki).name, rep_names{j})
            rep_idx(j) = ki;
            break;
        end
    end
    if rep_idx(j) == 0
        error('Representative orbit ''%s'' not found in benchmark_suite.mat', rep_names{j});
    end
end

if ~exist('results/sensitivity', 'dir')
    mkdir('results/sensitivity');
end

for j = 1:numel(rep_idx)
    orbit_idx = rep_idx(j);
    ref = load(fullfile('ref_trajectories', ...
               sprintf('ref_%s.mat', suite(orbit_idx).name)));

    run_idx      = 0;
    n_runs       = numel(delta_vec) * numel(N_vec) * numel(alpha_vec5);
    sens_results = repmat(struct('delta',[], 'N',[], 'alpha',[], 'gbl',[]), n_runs, 1);

    for delta_idx = 1:numel(delta_vec)
        for N_idx = 1:numel(N_vec)
            for alpha_idx = 1:numel(alpha_vec5)
                run_idx = run_idx + 1;
                fprintf('[%s/J6dragSRPmoon] delta=%2d N=%2d alpha=%.2f  (%d/%d) ... ', ...
                        suite(orbit_idx).name, ...
                        delta_vec(delta_idx), N_vec(N_idx), alpha_vec5(alpha_idx), ...
                        run_idx, n_runs);
                [gbl, ~, ~] = run_sweep_case(suite(orbit_idx), force_model, ...
                                              alpha_vec5(alpha_idx), ...
                                              delta_vec(delta_idx), N_vec(N_idx), ...
                                              ref.t_ref, ref.X_ref, C);
                sens_results(run_idx).delta = delta_vec(delta_idx);
                sens_results(run_idx).N     = N_vec(N_idx);
                sens_results(run_idx).alpha = alpha_vec5(alpha_idx);
                sens_results(run_idx).gbl   = gbl;
                fprintf('r_err=%.3e km\n', gbl.r_err_final);
            end
        end
    end

    orbit = suite(orbit_idx);
    fname = fullfile('results/sensitivity', ...
                     sprintf('sens_%s_%s.mat', orbit.name, force_model));
    save(fname, 'sens_results', 'delta_vec', 'N_vec', 'alpha_vec5', ...
         'orbit', 'force_model');
    fprintf('  -> Saved %s\n', fname);
end
