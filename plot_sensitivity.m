% plot_sensitivity.m
% Produces 3D scatter figures and tables from the sensitivity sweep results.
%
% Figures (no subplots):
%   One 3D scatter per representative orbit x FM (8 total)
%   Axes: alpha, N, delta  |  Color: log10(r_err_final)
%
% Tables (results/tables/):
%   sens_{name}_{fm}.tex — r_err_final for all (delta, N, alpha) combinations

%% Directories
sens_dir = fullfile('results', 'sensitivity');
fig_dir  = fullfile('results', 'figures', 'sensitivity');
tbl_dir  = fullfile('results', 'tables');
for d = {fig_dir, tbl_dir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

%% Setup
rep_orbits  = {'LEO-1', 'MEO-1', 'HEO-E07', 'GEO-1'};
force_model = 'J6dragSRPmoon';
delta_vec    = [4, 8, 16, 32];
N_vec        = [8, 16, 24, 32];
alpha_vec5   = [0, 0.25, 0.5, 0.75, 1.0];
n_delta      = numel(delta_vec);
n_N          = numel(N_vec);
n_alpha5     = numel(alpha_vec5);
[~, col5]    = get_alpha_colors();

% Shared marker encoding (same as plot_full)
mkr_shapes  = {'o','s','^','d'};       % one per delta value [4 8 16 32]
delta_enc   = [4, 8, 16, 32];
N_enc_all   = [4, 8, 16, 24, 32];
sz_enc_all  = [15, 25, 40, 60, 80];

%% Main loop
for oi = 1:numel(rep_orbits)
        name  = rep_orbits{oi};
        fm    = force_model;
        fpath = fullfile(sens_dir, sprintf('sens_%s_%s.mat', name, fm));
        if ~exist(fpath, 'file')
            warning('Missing: %s', fpath);
            continue
        end
        d = load(fpath);

        % Extract flat vectors directly from struct array
        r_err_flat = arrayfun(@(s) s.gbl.r_err_final,   d.sens_results);
        wt_flat    = arrayfun(@(s) s.gbl.wall_time,    d.sens_results);
        feval_flat = arrayfun(@(s) s.gbl.total_fevals, d.sens_results);
        alpha_flat = arrayfun(@(s) s.alpha,             d.sens_results);
        N_flat     = arrayfun(@(s) s.N,                 d.sens_results);
        delta_flat = arrayfun(@(s) s.delta,             d.sens_results);

        ok  = isfinite(r_err_flat);
        bad = ~ok;

        % ---- Unified 2D scatter figures (fevals vs r_err, wall_time vs r_err) --
        scatter_cfgs = {
            feval_flat, 'Total function evaluations', sprintf('%s_%s_eff',    name, fm);
            wt_flat,    'Wall time (s)',               sprintf('%s_%s_timing', name, fm)
        };
        for sc = 1:size(scatter_cfgs, 1)
            x_data  = scatter_cfgs{sc, 1};
            x_label = scatter_cfgs{sc, 2};
            sbase   = fullfile(fig_dir, scatter_cfgs{sc, 3});

            ok_s  = isfinite(x_data) & x_data > 0 & isfinite(r_err_flat) & r_err_flat > 0;
            bad_s = ~ok_s;

            fig_s = figure('Position', [50 50 900 550]);
            ax_s  = gca;
            hold(ax_s, 'on');
            set(ax_s, 'XScale', 'log', 'YScale', 'log');

            for ai = 1:n_alpha5
                has_any = false;
                for di = 1:n_delta
                    for ni = 1:n_N
                        mask = ok_s & abs(alpha_flat - alpha_vec5(ai)) < 0.01 ...
                                   & delta_flat == delta_vec(di) ...
                                   & N_flat     == N_vec(ni);
                        if ~any(mask), continue; end
                        sz = sz_enc_all(N_enc_all == N_vec(ni));
                        scatter(ax_s, x_data(mask), r_err_flat(mask), sz, col5(ai,:), ...
                                mkr_shapes{delta_enc == delta_vec(di)}, 'filled', ...
                                'HandleVisibility', 'off');
                        has_any = true;
                    end
                end
                if has_any
                    scatter(ax_s, NaN, NaN, 30, col5(ai,:), 'o', 'filled', ...
                            'DisplayName', sprintf('\\alpha=%.2g', alpha_vec5(ai)));
                end
            end

            % Shape/size key entries
            for di = 1:n_delta
                scatter(ax_s, NaN, NaN, 30, [0.4 0.4 0.4], ...
                        mkr_shapes{delta_enc == delta_vec(di)}, 'filled', ...
                        'DisplayName', sprintf('\\delta=%d', delta_vec(di)));
            end
            for ni = 1:n_N
                sz = sz_enc_all(N_enc_all == N_vec(ni));
                scatter(ax_s, NaN, NaN, sz, [0.4 0.4 0.4], 'o', 'filled', ...
                        'DisplayName', sprintf('N=%d', N_vec(ni)));
            end

            if any(bad_s)
                scatter(ax_s, NaN, NaN, 15, [0.6 0.6 0.6], 'x', ...
                        'DisplayName', 'failed');
            end

            xlabel(ax_s, x_label);
            ylabel(ax_s, 'Position error at t_{end} (km)');
            title(ax_s, sprintf('%s / %s', name, strrep(fm, '_', '\\_')));
            grid(ax_s, 'on'); box(ax_s, 'on');
            legend(ax_s, 'show', 'Location', 'northeast', 'FontSize', 7);
            save_fig(fig_s, sbase, true);
            fprintf('Saved scatter: %s\n', sbase);
        end

        % ---- Table -------------------------------------------------------
        % Reshape flat array to [n_delta x n_N x n_alpha5]
        % Loop order in sweep_sensitivity.m: delta outer -> N -> alpha inner
        err_cube = permute(reshape(r_err_flat, [n_alpha5, n_N, n_delta]), [3, 2, 1]);

        fname = fullfile(tbl_dir, sprintf('sens_%s_%s.tex', name, fm));
        fid   = fopen(fname, 'w');
        if fid < 0
            warning('Cannot write %s', fname);
            continue
        end
        fprintf(fid, '%% Sensitivity: %s / %s  (r_err_final in km)\n\n', name, fm);
        for ai = 1:n_alpha5
            fprintf(fid, '\\subsection*{$\\alpha = %.2g$}\n', alpha_vec5(ai));
            fprintf(fid, '\\begin{tabular}{c|%s}\n', repmat('r', 1, n_N));
            fprintf(fid, '$\\delta \\setminus N$');
            for ni = 1:n_N
                fprintf(fid, ' & %d', N_vec(ni));
            end
            fprintf(fid, ' \\\\\n\\hline\n');
            for di = 1:n_delta
                fprintf(fid, '%d', delta_vec(di));
                for ni = 1:n_N
                    val = err_cube(di, ni, ai);
                    if isfinite(val) && val > 0
                        fprintf(fid, ' & \\num{%.3e}', val);
                    else
                        fprintf(fid, ' & ---');
                    end
                end
                fprintf(fid, ' \\\\\n');
            end
            fprintf(fid, '\\end{tabular}\n\n');
        end
        fclose(fid);
        fprintf('Wrote r_err table: sens_%s_%s.tex\n', name, fm);

        % ---- wall_time and fevals tables ------------------------------------
        wt_cube    = permute(reshape(wt_flat,    [n_alpha5, n_N, n_delta]), [3, 2, 1]);
        feval_cube = permute(reshape(feval_flat, [n_alpha5, n_N, n_delta]), [3, 2, 1]);
        tbl_metrics = {wt_cube,    'wall time (s)',  'wtime', '%.2f'; ...
                       feval_cube, 'f-evals',        'fevals', '%.0f'};
        for tmi = 1:size(tbl_metrics,1)
            cube   = tbl_metrics{tmi,1};
            tlabel = tbl_metrics{tmi,2};
            tsuf   = tbl_metrics{tmi,3};
            tfmt   = tbl_metrics{tmi,4};
            tname  = fullfile(tbl_dir, sprintf('sens_%s_%s_%s.tex', name, fm, tsuf));
            fid2   = fopen(tname, 'w');
            if fid2 < 0, warning('Cannot write %s', tname); continue; end
            fprintf(fid2, '%% Sensitivity: %s / %s  (%s)\n\n', name, fm, tlabel);
            for ai = 1:n_alpha5
                fprintf(fid2, '\\subsection*{$\\alpha = %.2g$}\n', alpha_vec5(ai));
                fprintf(fid2, '\\begin{tabular}{c|%s}\n', repmat('r', 1, n_N));
                fprintf(fid2, '$\\delta \\setminus N$');
                for ni = 1:n_N
                    fprintf(fid2, ' & %d', N_vec(ni));
                end
                fprintf(fid2, ' \\\\\n\\hline\n');
                for di = 1:n_delta
                    fprintf(fid2, '%d', delta_vec(di));
                    for ni = 1:n_N
                        val = cube(di, ni, ai);
                        if isfinite(val) && val > 0
                            fprintf(fid2, [' & \\num{' tfmt '}'], val);
                        else
                            fprintf(fid2, ' & ---');
                        end
                    end
                    fprintf(fid2, ' \\\\\n');
                end
                fprintf(fid2, '\\end{tabular}\n\n');
            end
            fclose(fid2);
            fprintf('Wrote %s table: sens_%s_%s_%s.tex\n', tlabel, name, fm, tsuf);
        end
end

fprintf('\nplot_sensitivity done.\n');

% =========================================================================
function save_fig(fig, base, do_close)
figure(fig);
saveas(fig, [base '.png']);
try
    pos = get(gcf, 'Position');
    JGCD_W_PX = 2100;  % 7.0 in x 300 PPI (AIAA/JGCD full text width)
    JGCD_H_PX = round(JGCD_W_PX * pos(4) / pos(3));
    cleanfigure('targetResolution', [JGCD_W_PX, JGCD_H_PX]);
    matlab2tikz([base '.tikz'], 'width', '\linewidth', 'maxChunkLength', 500, 'showInfo', false);
catch e
    warning('tikz export failed for %s: %s', base, e.message);
end
if do_close
    close(fig);
end
end
