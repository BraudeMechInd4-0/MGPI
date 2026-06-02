% plot_alpha_sweep.m
% Produces all figures and tables from the alpha sweep results.
%
% Per-orbit figures (5 per orbit × 30 orbits = 150 total):
%   1. Per-node position error vs time                    [error_time/]
%   2. Per-node velocity error vs time                    [verr_time/]
%   3. Per-segment position error vs E (rad, unwrapped)   [error_E/]
%   4. d(err_pos)/dE vs E                                 [derrdE/]
%   5. Per-segment fevals vs E                            [fevals_E/]
%
% Global summary figures (8 runtime + 8 fevals + 4 error + 1 legend = 21 total):
%   6–9.   Wall time vs {e, h_perigee, inc, orbit name}   [runtime/]
%   10–13. Total fevals vs {e, h_perigee, inc, orbit name} [runtime/]
%   14–17. r_err_final vs {e, h_perigee, inc, orbit name} [runtime/]
%   18.    Shared legend                                   [root]
%
% Tables (results/tables/):
%   alpha_r_err_J6dragSRPmoon.tex
%   alpha_wtime_J6dragSRPmoon.tex
%   alpha_fevals_J6dragSRPmoon.tex

%% Parameters
t_zoom_start     = 0;     % days — left edge of error-vs-time window
t_zoom_end       = Inf;   % days — Inf = full span
close_after_save = true; % true: close each figure after saving
% orbit_names_filter: set before calling this script to restrict which orbits are plotted.
% Leave undefined or empty for all 30 orbits (default).
if ~exist('orbit_names_filter', 'var') || isempty(orbit_names_filter)
    orbit_names_filter = {};
end

[~, mu, ~, ~, ~, ~, ~, ~] = getgravc(84);

%% Directories
sweep_dir    = fullfile('results', 'sweep');
fig_et_dir   = fullfile('results', 'figures', 'alpha_sweep', 'error_time');
fig_vt_dir   = fullfile('results', 'figures', 'alpha_sweep', 'verr_time');
fig_eE_dir   = fullfile('results', 'figures', 'alpha_sweep', 'error_E');
fig_drdE_dir = fullfile('results', 'figures', 'alpha_sweep', 'derrdE');
fig_fE_dir   = fullfile('results', 'figures', 'alpha_sweep', 'fevals_E');
fig_rt_dir   = fullfile('results', 'figures', 'alpha_sweep', 'runtime');
fig_leg_dir  = fullfile('results', 'figures', 'alpha_sweep');
tbl_dir      = fullfile('results', 'tables');
for d = {fig_et_dir, fig_vt_dir, fig_eE_dir, fig_drdE_dir, fig_fE_dir, fig_rt_dir, tbl_dir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

%% Constants
orbit_names = {'LEO-1','LEO-2','MEO-1','GEO-1', ...
               'HEO-E01','HEO-E02','HEO-E03','HEO-E04','HEO-E05', ...
               'HEO-E06','HEO-E07','HEO-E08','HEO-E09','HEO-E10', ...
               'HEO-P01','HEO-P02','HEO-P03','HEO-P04','HEO-P05', ...
               'HEO-P06','HEO-P07','HEO-P08', ...
               'HEO-I01','HEO-I02','HEO-I03','HEO-I04','HEO-I05', ...
               'HEO-I06','HEO-I07','HEO-I08'};
n_orbits     = numel(orbit_names);
force_model  = 'J6dragSRPmoon';
alpha_vec    = 0 : 0.1 : 1.0;
n_alpha      = numel(alpha_vec);

% Colors: maxdistcolor palette, fixed seed — all 11 α values
col_all   = get_alpha_colors();

% Runtime/fevals global figures: all 11 α values
sub_idx   = 1:n_alpha;
alpha_sub = alpha_vec;
col_sub   = col_all;

%% Table storage [n_orbits x n_alpha]
r_err_tbl  = NaN(n_orbits, n_alpha);
wtime_tbl  = NaN(n_orbits, n_alpha);
fevals_tbl = NaN(n_orbits, n_alpha);

%% Orbit metadata (loaded from sweep files — orbit struct embedded in each .mat)
suite_e        = NaN(n_orbits, 1);
suite_h_peri   = NaN(n_orbits, 1);
suite_inc      = NaN(n_orbits, 1);
suite_sweep    = cell(n_orbits, 1);

%% Main loop: 5 figures per orbit
for oi = 1:n_orbits
    name = orbit_names{oi};
    if ~isempty(orbit_names_filter) && ~any(strcmp(orbit_names_filter, name))
        continue
    end
    fpath = fullfile(sweep_dir, sprintf('sweep_%s_%s.mat', name, force_model));
    if ~exist(fpath, 'file')
        warning('Missing: %s', fpath);
        continue
    end
    d = load(fpath);

    % Collect metadata from embedded orbit struct
    if isfield(d, 'orbit')
        suite_e(oi)      = d.orbit.e;
        suite_h_peri(oi) = d.orbit.h_perigee;
        suite_inc(oi)    = d.orbit.inc;
        suite_sweep{oi}  = d.orbit.sweep;
    end

    for a = 1:n_alpha
        if d.global_res(a).error_code ~= 0, continue; end
        r_err_tbl(oi, a)  = d.global_res(a).r_err_final;
        wtime_tbl(oi, a)  = d.global_res(a).wall_time;
        fevals_tbl(oi, a) = d.global_res(a).total_fevals;
    end

    fig_title = sprintf('%s / %s', name, strrep(force_model, '_', '\\_'));

    % Figure 1: per-node position error vs time
    fig = figure('Position', [50 50 1200 400]);
    ax  = gca;
    hold(ax, 'on');
    set(ax, 'YScale', 'log');
    t_max_days = 0;
    for a = 1:n_alpha
        if d.global_res(a).error_code ~= 0, continue; end
        t_days = d.nod_res{a}.t / 86400;
        plot(ax, t_days, d.nod_res{a}.r_err, ...
             'Color', col_all(a,:), 'LineWidth', 0.7, ...
             'DisplayName', sprintf('\\alpha=%.1f', alpha_vec(a)));
        t_max_days = max(t_max_days, t_days(end));
    end
    xlim(ax, [t_zoom_start, min(t_zoom_end, t_max_days)]);
    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Position error (km)');
    title(ax, fig_title);
    grid(ax, 'on'); box(ax, 'on');
    base = fullfile(fig_et_dir, sprintf('%s_%s', name, force_model));
    save_fig(fig, base, '\linewidth', close_after_save);

    % Figure 2: per-node velocity error vs time
    fig = figure('Position', [50 50 1200 400]);
    ax  = gca;
    hold(ax, 'on');
    set(ax, 'YScale', 'log');
    t_max_days = 0;
    for a = 1:n_alpha
        if d.global_res(a).error_code ~= 0, continue; end
        t_days = d.nod_res{a}.t / 86400;
        plot(ax, t_days, d.nod_res{a}.v_err, ...
             'Color', col_all(a,:), 'LineWidth', 0.7, ...
             'DisplayName', sprintf('\\alpha=%.1f', alpha_vec(a)));
        t_max_days = max(t_max_days, t_days(end));
    end
    xlim(ax, [t_zoom_start, min(t_zoom_end, t_max_days)]);
    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Velocity error (km/s)');
    title(ax, fig_title);
    grid(ax, 'on'); box(ax, 'on');
    base = fullfile(fig_vt_dir, sprintf('%s_%s', name, force_model));
    save_fig(fig, base, '\linewidth', close_after_save);

    % Figure 3: per-segment position error vs E (unwrapped)
    fig = figure('Position', [50 50 1200 400]);
    ax  = gca;
    hold(ax, 'on');
    set(ax, 'YScale', 'log');
    for a = 1:n_alpha
        if d.global_res(a).error_code ~= 0, continue; end
        E_raw  = d.seg_res{a}.E_osc;
        y_data = d.seg_res{a}.r_err;
        valid  = ~isnan(E_raw);
        if ~any(valid), continue; end
        E_uw = NaN(size(E_raw));
        E_uw(valid) = unwrap(E_raw(valid));
        plot(ax, E_uw, y_data, ...
             'Color', col_all(a,:), 'LineWidth', 0.7, ...
             'DisplayName', sprintf('\\alpha=%.1f', alpha_vec(a)));
    end
    xlabel(ax, 'E (rad)');
    ylabel(ax, 'Segment position error (km)');
    title(ax, fig_title);
    grid(ax, 'on'); box(ax, 'on');
    base = fullfile(fig_eE_dir, sprintf('%s_%s', name, force_model));
    save_fig(fig, base, '\linewidth', close_after_save);

    % Figure 4: d(err_pos)/dE vs E (unwrapped)
    fig = figure('Position', [50 50 1200 400]);
    ax  = gca;
    hold(ax, 'on');
    set(ax, 'YScale', 'log');
    for a = 1:n_alpha
        if d.global_res(a).error_code ~= 0, continue; end
        E_raw  = d.seg_res{a}.E_osc;
        e_s    = d.seg_res{a}.e_osc;
        r_s    = d.seg_res{a}.r_start;
        v_s    = d.seg_res{a}.v_start;
        a_osc  = 1 ./ (2./r_s - v_s.^2 ./ mu);
        y_data = d.seg_res{a}.v_err .* (1 - e_s .* cos(E_raw)) .* sqrt(a_osc.^3 ./ mu);
        valid  = ~isnan(E_raw);
        if ~any(valid), continue; end
        E_uw = NaN(size(E_raw));
        E_uw(valid) = unwrap(E_raw(valid));
        plot(ax, E_uw, y_data, ...
             'Color', col_all(a,:), 'LineWidth', 0.7, ...
             'DisplayName', sprintf('\\alpha=%.1f', alpha_vec(a)));
    end
    xlabel(ax, 'E (rad)');
    ylabel(ax, 'd(r_{err})/dE (km/rad)');
    title(ax, fig_title);
    grid(ax, 'on'); box(ax, 'on');
    base = fullfile(fig_drdE_dir, sprintf('%s_%s', name, force_model));
    save_fig(fig, base, '\linewidth', close_after_save);

    % Figure 5: per-segment fevals vs E (unwrapped)
    fig = figure('Position', [50 50 1200 400]);
    ax  = gca;
    hold(ax, 'on');
    for a = 1:n_alpha
        if d.global_res(a).error_code ~= 0, continue; end
        E_raw  = d.seg_res{a}.E_osc;
        y_data = d.seg_res{a}.fevals;
        valid  = ~isnan(E_raw);
        if ~any(valid), continue; end
        E_uw = NaN(size(E_raw));
        E_uw(valid) = unwrap(E_raw(valid));
        plot(ax, E_uw, y_data, ...
             'Color', col_all(a,:), 'LineWidth', 0.7, ...
             'DisplayName', sprintf('\\alpha=%.1f', alpha_vec(a)));
    end
    xlabel(ax, 'E (rad)');
    ylabel(ax, 'Fevals per segment');
    title(ax, fig_title);
    grid(ax, 'on'); box(ax, 'on');
    base = fullfile(fig_fE_dir, sprintf('%s_%s', name, force_model));
    save_fig(fig, base, '\linewidth', close_after_save);

    fprintf('Saved per-orbit figures: %s\n', name);
end

%% Global runtime/fevals figures
% Only 5 α values to keep runtime plots readable
% Each metric (wall time, fevals) gets 4 figures: ecc, perigee, inc, anchors

make_runtime_figure(suite_sweep, suite_e,      wtime_tbl,  sub_idx, col_sub, alpha_sub, ...
    'ecc', 'eccentricity', 'Total wall time (s)', 'wtime_ecc',   fig_rt_dir, close_after_save, force_model);
make_runtime_figure(suite_sweep, suite_h_peri, wtime_tbl,  sub_idx, col_sub, alpha_sub, ...
    'perigee', 'Perigee altitude (km)', 'Total wall time (s)', 'wtime_perigee', fig_rt_dir, close_after_save, force_model);
make_runtime_figure(suite_sweep, suite_inc,    wtime_tbl,  sub_idx, col_sub, alpha_sub, ...
    'inc', 'Inclination (deg)', 'Total wall time (s)', 'wtime_inc',   fig_rt_dir, close_after_save, force_model);
make_runtime_figure_anchors(suite_sweep, orbit_names, wtime_tbl,  sub_idx, col_sub, alpha_sub, ...
    'Total wall time (s)', 'wtime_anchors',   fig_rt_dir, close_after_save, force_model);

make_runtime_figure(suite_sweep, suite_e,      fevals_tbl, sub_idx, col_sub, alpha_sub, ...
    'ecc', 'eccentricity', 'Total fevals', 'fevals_ecc',   fig_rt_dir, close_after_save, force_model);
make_runtime_figure(suite_sweep, suite_h_peri, fevals_tbl, sub_idx, col_sub, alpha_sub, ...
    'perigee', 'Perigee altitude (km)', 'Total fevals', 'fevals_perigee', fig_rt_dir, close_after_save, force_model);
make_runtime_figure(suite_sweep, suite_inc,    fevals_tbl, sub_idx, col_sub, alpha_sub, ...
    'inc', 'Inclination (deg)', 'Total fevals', 'fevals_inc',   fig_rt_dir, close_after_save, force_model);
make_runtime_figure_anchors(suite_sweep, orbit_names, fevals_tbl, sub_idx, col_sub, alpha_sub, ...
    'Total fevals', 'fevals_anchors',   fig_rt_dir, close_after_save, force_model);

make_runtime_figure(suite_sweep, suite_e,      r_err_tbl,  sub_idx, col_sub, alpha_sub, ...
    'ecc', 'eccentricity', 'Position error at t_{end} (km)', 'rerr_ecc',   fig_rt_dir, close_after_save, force_model, true);
make_runtime_figure(suite_sweep, suite_h_peri, r_err_tbl,  sub_idx, col_sub, alpha_sub, ...
    'perigee', 'Perigee altitude (km)', 'Position error at t_{end} (km)', 'rerr_perigee', fig_rt_dir, close_after_save, force_model, true);
make_runtime_figure(suite_sweep, suite_inc,    r_err_tbl,  sub_idx, col_sub, alpha_sub, ...
    'inc', 'Inclination (deg)', 'Position error at t_{end} (km)', 'rerr_inc',   fig_rt_dir, close_after_save, force_model, true);
make_runtime_figure_anchors(suite_sweep, orbit_names, r_err_tbl,  sub_idx, col_sub, alpha_sub, ...
    'Position error at t_{end} (km)', 'rerr_anchors', fig_rt_dir, close_after_save, force_model, true);

%% Shared legend (all 11 α, maxdistcolor)
fig_leg = figure('Position', [50 50 1400 100]);
ax_leg  = axes(fig_leg);
hold(ax_leg, 'on');
for a = 1:n_alpha
    plot(ax_leg, NaN, NaN, 'Color', col_all(a,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('\\alpha = %.1f', alpha_vec(a)));
end
legend(ax_leg, 'show', 'Orientation', 'horizontal', 'Location', 'best');
axis(ax_leg, 'off');
base_leg = fullfile(fig_leg_dir, 'legend');
save_fig(fig_leg, base_leg, [], close_after_save);
fprintf('Saved legend\n');

%% Tables
write_alpha_table(r_err_tbl,  orbit_names, alpha_vec, 'r\_err (km)',   force_model, tbl_dir, ...
    sprintf('alpha_r_err_%s.tex',  force_model), false);
write_alpha_table(wtime_tbl,  orbit_names, alpha_vec, 'wall time (s)', force_model, tbl_dir, ...
    sprintf('alpha_wtime_%s.tex',  force_model), false);
write_alpha_table(fevals_tbl, orbit_names, alpha_vec, 'f-evals',       force_model, tbl_dir, ...
    sprintf('alpha_fevals_%s.tex', force_model), true);

fprintf('\nplot_alpha_sweep done.\n');

% =========================================================================
function make_runtime_figure(sweep_cell, x_data, y_tbl, sub_idx, col, alpha_sub, ...
                              sweep_tag, x_label, y_label, fname_base, out_dir, do_close, fm, log_y)
% Plot one global figure: y_tbl(:, sub_idx) vs x_data for orbits whose
% sweep field matches sweep_tag. Also includes molniya in ecc plot.
if nargin < 14, log_y = false; end
mask = strcmp(sweep_cell, sweep_tag);
if strcmp(sweep_tag, 'ecc')
    mask = mask | strcmp(sweep_cell, 'molniya');
end
if ~any(mask), return; end

x = x_data(mask);
y = y_tbl(mask, sub_idx);          % (n_filtered × n_sub)
[x_sorted, ix] = sort(x);
y_sorted = y(ix, :);

fig = figure('Position', [50 50 900 420]);
ax  = gca;
hold(ax, 'on');
if log_y, set(ax, 'YScale', 'log'); end
for j = 1:numel(alpha_sub)
    plot(ax, x_sorted, y_sorted(:, j), ...
         'Color', col(j,:), 'LineWidth', 1.2, ...
         'Marker', 'o', 'MarkerSize', 5, ...
         'DisplayName', sprintf('\\alpha=%.2f', alpha_sub(j)));
end
xlabel(ax, x_label);
ylabel(ax, y_label);
title(ax, sprintf('%s vs %s — %s', y_label, sweep_tag, strrep(fm, '_', '\\_')));
grid(ax, 'on'); box(ax, 'on');
legend(ax, 'show', 'Location', 'best');
base = fullfile(out_dir, [fname_base '_' fm]);
save_fig(fig, base, '\linewidth', do_close);
fprintf('Saved runtime figure: %s\n', fname_base);
end

% =========================================================================
function make_runtime_figure_anchors(sweep_cell, orbit_names, y_tbl, sub_idx, col, alpha_sub, ...
                                      y_label, fname_base, out_dir, do_close, fm, log_y)
if nargin < 12, log_y = false; end
mask = strcmp(sweep_cell, 'anchor');
if ~any(mask), return; end

names_anch = orbit_names(mask);
y = y_tbl(mask, sub_idx);          % (n_anchors × n_sub)

n_anch  = numel(names_anch);
x_pos   = (0:n_anch-1) * 3 + 1;   % wider spacing between orbit groups

fig = figure('Position', [50 50 700 420]);
ax  = gca;
hold(ax, 'on');
if log_y, set(ax, 'YScale', 'log'); end
for j = 1:numel(alpha_sub)
    bar_offset = (j - (numel(alpha_sub)+1)/2) * 0.18;   % per-alpha spacing within group
    scatter(ax, x_pos + bar_offset, y(:, j), 40, col(j,:), 'filled', ...
            'DisplayName', sprintf('\\alpha=%.2f', alpha_sub(j)));
end
set(ax, 'XTick', x_pos, 'XTickLabel', names_anch);
ylabel(ax, y_label);
title(ax, sprintf('%s — circular anchors — %s', y_label, strrep(fm, '_', '\\_')));
grid(ax, 'on'); box(ax, 'on');
legend(ax, 'show', 'Location', 'best');
base = fullfile(out_dir, [fname_base '_' fm]);
save_fig(fig, base, '\linewidth', do_close);
fprintf('Saved anchor runtime figure: %s\n', fname_base);
end

% =========================================================================
function save_fig(fig, base, tikz_width, do_close)
figure(fig);
saveas(fig, [base '.png']);
try
    pos = get(gcf, 'Position');
    JGCD_W_PX = 2100;  % 7.0 in x 300 PPI (AIAA/JGCD full text width)
    JGCD_H_PX = round(JGCD_W_PX * pos(4) / pos(3));
    cleanfigure('targetResolution', [JGCD_W_PX, JGCD_H_PX]);
    if isempty(tikz_width)
        matlab2tikz([base '.tikz'], 'width', '\linewidth', 'maxChunkLength', 500, 'showInfo', false);
    else
        matlab2tikz([base '.tikz'], 'width', tikz_width, 'maxChunkLength', 500, 'showInfo', false);
    end
catch e
    warning('tikz export failed for %s: %s', base, e.message);
end
if do_close
    close(fig);
end
end

% =========================================================================
function write_alpha_table(data, orbit_names, alpha_vec, metric_label, fm_label, tbl_dir, fname, is_integer)
n_orbits = size(data, 1);
n_alpha  = size(data, 2);

fprintf('\n=== %s / %s ===\n', strrep(metric_label, '\', ''), fm_label);
fprintf('%-10s', 'Orbit');
for a = 1:n_alpha
    fprintf('  %9.2f', alpha_vec(a));
end
fprintf('\n');
for oi = 1:n_orbits
    fprintf('%-10s', orbit_names{oi});
    for a = 1:n_alpha
        if isnan(data(oi, a))
            fprintf('  %9s', '---');
        elseif is_integer
            fprintf('  %9d', round(data(oi, a)));
        else
            fprintf('  %9.3e', data(oi, a));
        end
    end
    fprintf('\n');
end

fid = fopen(fullfile(tbl_dir, fname), 'w');
if fid < 0
    warning('Cannot write %s', fname);
    return
end
fprintf(fid, '\\begin{tabular}{l|%s}\n\\hline\n', repmat('r', 1, n_alpha));
fprintf(fid, 'Orbit');
for a = 1:n_alpha
    fprintf(fid, ' & $\\alpha{=}%.1f$', alpha_vec(a));
end
fprintf(fid, ' \\\\\n\\hline\n');
for oi = 1:n_orbits
    fprintf(fid, '%s', orbit_names{oi});
    for a = 1:n_alpha
        if isnan(data(oi, a))
            fprintf(fid, ' & ---');
        elseif is_integer
            fprintf(fid, ' & %d', round(data(oi, a)));
        else
            fprintf(fid, ' & \\num{%.3e}', data(oi, a));
        end
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\hline\n\\end{tabular}\n');
fclose(fid);
fprintf('Wrote %s\n', fullfile(tbl_dir, fname));
end
