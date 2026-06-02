% run_all_sweeps.m
%
% Full pipeline: initialises the benchmark suite then runs all three sweeps.
%
%   Step 0   init_benchmark_suite   30 orbit structs + 30 reference trajectories
%                                   (ode78 at 1e-13, 30 days each)
%                                   ~30-120 min depending on hardware
%
%   Stage 1  sweep_alpha            alpha sweep, delta=16 N=16, 30 days
%                                   30 orbits x 11 alpha = 330 MPCM runs
%
%   Stage 2  sweep_sensitivity      delta x N x alpha sensitivity, 30 days
%                                   4 rep orbits x 4delta x 4N x 5alpha = 320 runs
%
%   Stage 3  sweep_full             timing + accuracy, 7-day slice
%                                   30 orbits x (4delta x 4N x 11alpha) = 5280 MPCM runs
%                                   + 30 x 4 ode45 baseline runs
%
%   Stage 4  plot_alpha_sweep       figures + tables from sweep_alpha results
%            plot_sensitivity       figures + tables from sweep_sensitivity results
%            plot_full              figures from sweep_full results
%
% Estimated total runtime: 1-4 days (hardware-dependent; HEO orbits are slow).
% Results are saved incrementally — each stage is safe to resume independently.
% Run test_sweep_alpha.m first to verify the setup on LEO-1 alone.

fprintf('=== Pipeline start: %s ===\n\n', char(datetime('now')));
t_pipeline = tic;

%% Step 0: generate benchmark suite and reference trajectories
n_ref_expected = 30;
ref_files = dir(fullfile('ref_trajectories', 'ref_*.mat'));
if exist('benchmark_suite.mat', 'file') == 2 && ...
        exist('ref_trajectories', 'dir') == 7 && ...
        numel(ref_files) >= n_ref_expected
    fprintf('Step 0: benchmark_suite.mat and %d ref files found — skipping init.\n\n', ...
        numel(ref_files));
    t0_min = 0;
else
    fprintf('--- Step 0 / 3: init_benchmark_suite ---\n');
    t0 = tic;
    init_benchmark_suite
    t0_min = toc(t0)/60;
    fprintf('\nStep 0 done in %.1f min\n\n', t0_min);
end

%% Stage 1: alpha sweep (sweep_alpha.m)
fprintf('--- Stage 1 / 4: sweep_alpha ---\n');
t1 = tic;
sweep_alpha
fprintf('\nStage 1 done in %.1f min\n\n', toc(t1)/60);

%% Stage 2: sensitivity sweep (sweep_sensitivity.m)
fprintf('--- Stage 2 / 4: sweep_sensitivity ---\n');
t2 = tic;
sweep_sensitivity
fprintf('\nStage 2 done in %.1f min\n\n', toc(t2)/60);

%% Stage 3: full timing sweep (sweep_full.m)
fprintf('--- Stage 3 / 4: sweep_full ---\n');
t3 = tic;
sweep_full
fprintf('\nStage 3 done in %.1f min\n\n', toc(t3)/60);

%% Stage 4: plots and tables
fprintf('--- Stage 4 / 4: plots ---\n');
t4 = tic;
plot_alpha_sweep
plot_sensitivity
plot_full
fprintf('\nStage 4 done in %.1f min\n\n', toc(t4)/60);

%% Summary
fprintf('=== Pipeline complete: %s ===\n', char(datetime('now')));
if t0_min > 0
    fprintf('  Step 0  (init_benchmark_suite): %.1f min\n', t0_min);
end
fprintf('  Stage 1 (sweep_alpha):          %.1f min\n', toc(t1)/60);
fprintf('  Stage 2 (sweep_sensitivity):    %.1f min\n', toc(t2)/60);
fprintf('  Stage 3 (sweep_full):           %.1f min\n', toc(t3)/60);
fprintf('  Stage 4 (plots):                %.1f min\n', toc(t4)/60);
fprintf('  Total:                          %.2f h\n',   toc(t_pipeline)/3600);
