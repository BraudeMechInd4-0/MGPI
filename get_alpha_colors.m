function [col11, col5] = get_alpha_colors()
% Returns fixed color matrices for all MGPI alpha sweep figures.
%
%   col11 : 11x3 — one color per alpha value 0.0, 0.1, ..., 1.0
%   col5  :  5x3 — colors for sensitivity alpha = {0, 0.25, 0.5, 0.75, 1.0}
%              alpha=0.0, 0.5, 1.0 match their col11 counterparts
%              alpha=0.25 -> black
%              alpha=0.75 -> yellow  (not in maxdistcolor — avoids light print issue)
%
% Colors are generated once with a fixed rng seed for reproducibility.

rng(42);
col11 = maxdistcolor(11, 'Lmin', 0.45, 'Lmax', 0.7, 'Cmin', 0.1, 'sort', 'hue');

col5 = [
    col11(1,  :);   % alpha = 0.00
    0    0    0;    % alpha = 0.25 — black
    col11(6,  :);   % alpha = 0.50
    1    1    0;    % alpha = 0.75 — yellow
    col11(11, :);   % alpha = 1.00
];
end