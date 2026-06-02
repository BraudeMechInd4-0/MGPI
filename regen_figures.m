% regen_figures.m
% Re-plots all figures at JGCD 300 PPI, checks data-point count, rebuilds processed tikz.
% Equivalent to Stage 4 of run_all_sweeps.m, plus verification and post-processing.

fprintf('--- Stage 4: plots ---\n'); t = tic;
plot_alpha_sweep
plot_sensitivity
plot_full
fprintf('Plots done in %.1f min\n\n', toc(t)/60);

%% Point-count and max-line-length check
sample = 'results/figures/alpha_sweep/error_time/LEO-1_J6dragSRPmoon.tikz';
content = fileread(sample);
n_pts = numel(strfind(content, '\\'));
fprintf('LEO-1 error_time data rows: %d  (target ~15,000 per figure)\n', n_pts);
fid = fopen(sample, 'r'); maxLen = 0;
while ~feof(fid), ln = fgetl(fid); if ischar(ln), maxLen = max(maxLen, length(ln)); end; end
fclose(fid);
fprintf('Longest line: %d chars  (must be < 200,000)\n\n', maxLen);

%% Rebuild processed tikz folder (PROCESS_TIKZ=True: patch + copy all tikz files)
fprintf('--- gen_figures_tex.py (PROCESS_TIKZ=True) ---\n');
[status, out] = system('python -c "import gen_figures_tex as g; g.PROCESS_TIKZ=True; g.main()"');
fprintf('%s\n', out);
if status ~= 0
    warning('gen_figures_tex.py exited with status %d', status);
end
fprintf('Done. Copy results/figures_processed/ -> Figures/ in LaTeX project.\n');
