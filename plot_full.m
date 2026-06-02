% plot_full.m
% Efficiency comparison figures from sweep_full results.
%
% Purpose: for each orbit, scatter (fevals, r_err_final) for all MPCM
% δ×N×α combinations vs ode45 at 4 tolerance settings.
% Answers: at a given accuracy target, which method needs fewer function evaluations?
%
% Figures (saved to results/figures/full/):
%   One scatter per orbit (30 total): fevals vs r_err_final
%     MPCM points colored by α (parula(11) gradient)
%     ode45 points marked with distinct symbols

%% Directories
full_dir = fullfile('results', 'full');
fig_dir  = fullfile('results', 'figures', 'full');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Constants
force_model = 'J6dragSRPmoon';
alpha_vec   = 0 : 0.1 : 1.0;
n_alpha     = numel(alpha_vec);
col         = get_alpha_colors();

% Shared marker encoding (same as plot_sensitivity)
mkr_shapes  = {'o','s','^','d'};           % one per delta value [4 8 16 32]
delta_enc   = [4, 8, 16, 32];             % must match sweep_full delta_vec
N_enc_all   = [4, 8, 16, 24, 32];         % all N values across both sweeps
sz_enc_all  = [15, 25, 40, 60, 80];       % corresponding marker sizes

orbit_names = {'LEO-1','LEO-2','MEO-1','GEO-1', ...
               'HEO-E01','HEO-E02','HEO-E03','HEO-E04','HEO-E05', ...
               'HEO-E06','HEO-E07','HEO-E08','HEO-E09','HEO-E10', ...
               'HEO-P01','HEO-P02','HEO-P03','HEO-P04','HEO-P05', ...
               'HEO-P06','HEO-P07','HEO-P08', ...
               'HEO-I01','HEO-I02','HEO-I03','HEO-I04','HEO-I05', ...
               'HEO-I06','HEO-I07','HEO-I08'};

ode45_markers = {'s','d','^','v'};   % one per tolerance setting
ode45_color   = [0.8 0 0];          % red for all ode45 points

close_after_save = true;

%% Main loop
for oi = 1:numel(orbit_names)
    name  = orbit_names{oi};
    fpath = fullfile(full_dir, sprintf('full_%s_%s.mat', name, force_model));
    if ~exist(fpath, 'file')
        warning('Missing: %s', fpath);
        continue
    end
    d = load(fpath);

    % Extract MPCM data
    alpha_flat  = arrayfun(@(s) s.alpha,             d.full_results);
    fevals_flat = arrayfun(@(s) s.gbl.total_fevals,  d.full_results);
    rerr_flat   = arrayfun(@(s) s.gbl.r_err_final,   d.full_results);
    ec_flat     = arrayfun(@(s) s.gbl.error_code,    d.full_results);
    delta_flat  = arrayfun(@(s) s.delta,             d.full_results);
    N_flat      = arrayfun(@(s) s.N,                 d.full_results);

    ok  = (ec_flat == 0) & isfinite(rerr_flat) & (rerr_flat > 0) & (fevals_flat > 0);
    bad = ~ok;

    delta_vals = unique(delta_flat);
    N_vals     = unique(N_flat);

    fig = figure('Position', [50 50 900 550]);
    ax  = gca;
    hold(ax, 'on');
    set(ax, 'XScale', 'log', 'YScale', 'log');

    % MPCM points — triple loop: color=alpha, shape=delta, size=N
    for ai = 1:n_alpha
        has_any = false;
        for di = 1:numel(delta_vals)
            for ni = 1:numel(N_vals)
                mask = ok & abs(alpha_flat - alpha_vec(ai)) < 0.01 ...
                          & delta_flat == delta_vals(di) ...
                          & N_flat     == N_vals(ni);
                if ~any(mask), continue; end
                sz = sz_enc_all(N_enc_all == N_vals(ni));
                scatter(ax, fevals_flat(mask), rerr_flat(mask), sz, col(ai,:), ...
                        mkr_shapes{delta_enc == delta_vals(di)}, 'filled', ...
                        'HandleVisibility', 'off');
                has_any = true;
            end
        end
        if has_any  % dummy point for alpha legend entry
            scatter(ax, NaN, NaN, 30, col(ai,:), 'o', 'filled', ...
                    'DisplayName', sprintf('\\alpha=%.1f', alpha_vec(ai)));
        end
    end

    % Shape/size key as additional legend entries
    for di = 1:numel(delta_vals)
        scatter(ax, NaN, NaN, 30, [0.4 0.4 0.4], ...
                mkr_shapes{delta_enc == delta_vals(di)}, 'filled', ...
                'DisplayName', sprintf('\\delta=%d', delta_vals(di)));
    end
    for ni = 1:numel(N_vals)
        sz = sz_enc_all(N_enc_all == N_vals(ni));
        scatter(ax, NaN, NaN, sz, [0.4 0.4 0.4], 'o', 'filled', ...
                'DisplayName', sprintf('N=%d', N_vals(ni)));
    end

    % Failed MPCM points
    if any(bad) && any(ok)
        max_feval = max(fevals_flat(ok), [], 'all');
        scatter(ax, ones(sum(bad),1)*max_feval*1.05, ones(sum(bad),1)*1e-15, ...
                15, [0.6 0.6 0.6], 'x', 'DisplayName', 'MPCM failed');
    end

    % ode45 points
    for ti = 1:numel(d.ode45_results)
        o45 = d.ode45_results(ti);
        if ~isfinite(o45.r_err_final) || o45.r_err_final <= 0 || o45.fevals <= 0
            continue
        end
        tol_lbl = sprintf('ode45 %.0e/%.0e', o45.tol_rel, o45.tol_abs);
        mk = ode45_markers{min(ti, numel(ode45_markers))};
        scatter(ax, o45.fevals, o45.r_err_final, 80, ode45_color, mk, ...
                'LineWidth', 2, 'DisplayName', tol_lbl);
    end

    xlabel(ax, 'Total function evaluations');
    ylabel(ax, 'Position error at t_{end} (km)');
    title(ax, sprintf('%s / %s — efficiency', name, strrep(force_model, '_', '\\_')));
    grid(ax, 'on'); box(ax, 'on');
    legend(ax, 'show', 'Location', 'northeast', 'FontSize', 7);

    base = fullfile(fig_dir, sprintf('full_%s_%s', name, force_model));
    save_fig(fig, base, close_after_save);
    fprintf('Saved: %s\n', base);
end

fprintf('\nplot_full done.\n');

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
