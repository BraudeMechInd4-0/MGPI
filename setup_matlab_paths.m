% setup_matlab_paths.m
%
% Edit the five path variables below to point to your local copies of each
% dependency, then run this script once. You can also call it from startup.m
% so the paths are set automatically each MATLAB session.
%
% See README.md for download locations for each dependency.

SGP4_PATH     = 'C:/path/to/fundamentals-of-astrodynamics/software/matlab/SGP4';
OPQ_ORTH_PATH = 'C:/path/to/OPQ/OPQ_orthpol';
OPQ_QUAD_PATH = 'C:/path/to/OPQ/OPQ_quadrature';
M2TIKZ_PATH   = 'C:/path/to/matlab2tikz';
TBL2TEX_PATH  = 'C:/path/to/table2latex';

addpath(SGP4_PATH);
addpath(OPQ_ORTH_PATH);
addpath(OPQ_QUAD_PATH);
addpath(M2TIKZ_PATH);
addpath(TBL2TEX_PATH);

fprintf('MGPI paths configured.\n');