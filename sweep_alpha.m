% sweep_alpha.m  —  Stage 5.1
%
% α sweep: α ∈ {0, 0.1, …, 1.0} at fixed δ=16, N=16
% across all 30 orbits × 1 force model (J6dragSRPmoon).
%
% Output: results/sweep/sweep_<name>_J6dragSRPmoon.mat per orbit
%   alpha_vec   : 11×1  swept α values
%   global_res  : 11×1 struct  (see run_sweep_case.m for field names)
%   seg_res     : 11×1 cell of per-segment structs
%   nod_res     : 11×1 cell of per-node structs
%   orbit       : suite element
%   force_model : 'J6dragSRPmoon'
%   delta, N_nodes : fixed parameters used

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

alpha_vec    = (0 : 0.1 : 1.0)';
force_model  = 'J6dragSRPmoon';
delta        = 16;
N_nodes      = 16;

if ~exist('results/sweep', 'dir')
    mkdir('results/sweep');
end

for orbit_idx = 1:numel(suite)
    ref = load(fullfile('ref_trajectories', ...
               sprintf('ref_%s.mat', suite(orbit_idx).name)));

    n_alpha    = numel(alpha_vec);
    global_res = repmat(struct('alpha',[], 'r_err_final',[], 'v_err_final',[], ...
                               'total_iters',[], 'total_fevals',[], ...
                               'wall_time',[], 'error_code',[]), n_alpha, 1);
    seg_res = cell(n_alpha, 1);
    nod_res = cell(n_alpha, 1);

    for alpha_idx = 1:n_alpha
        fprintf('[%s/J6dragSRPmoon] alpha=%.1f ... ', ...
                suite(orbit_idx).name, alpha_vec(alpha_idx));
        [global_res(alpha_idx), seg_res{alpha_idx}, nod_res{alpha_idx}] = ...
            run_sweep_case(suite(orbit_idx), force_model, alpha_vec(alpha_idx), ...
                           delta, N_nodes, ref.t_ref, ref.X_ref, C);
        fprintf('r_err=%.3e km  fevals=%d  wall=%.1fs\n', ...
                global_res(alpha_idx).r_err_final, ...
                global_res(alpha_idx).total_fevals, ...
                global_res(alpha_idx).wall_time);
    end

    orbit = suite(orbit_idx);
    fname = fullfile('results/sweep', sprintf('sweep_%s_%s.mat', orbit.name, force_model));
    save(fname, 'alpha_vec', 'global_res', 'seg_res', 'nod_res', ...
         'orbit', 'force_model', 'delta', 'N_nodes');
    fprintf('  -> Saved %s\n', fname);
end
