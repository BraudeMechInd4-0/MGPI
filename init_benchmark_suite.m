% init_benchmark_suite.m
%
% Initialises the MGPI orbital benchmark suite (v3 — 30 orbits, single FM).
%
% Creates:
%   benchmark_suite.mat     - suite struct array (30 elements)
%   ref_trajectories/       - 30 reference .mat files
%
% Each reference file: t_ref (Nx1 s), X_ref (Nx6 km/km/s), orbit, force_model
% Force model: 'J6dragSRPmoon' -> J2-J6 + drag + SRP + lunar 3rd body
% Integrator: ode78 at RelTol=1e-13, AbsTol=1e-13
% Time grid:  T_orbital/2048 per orbit
% Duration:   30 days

%% Constants
[~, mu, Re, ~, ~, ~, ~, ~] = getgravc(84);
J = [0, 1.08262617385222e-3, -2.53241051856772e-6, -1.61989759991697e-6, ...
        -2.27753590730836e-7, 5.40666576283813e-7];  % Vallado FADA

CD_drag  = 2.2;
A_drag   = 4e-6;       % km²
m_drag   = 260;        % kg
jdepoch  = 2461042.0;  % Jan 1, 2026 12:00 TT
rhoSRP   = 4.56e-6;    % N/m²  (Vallado, solar pressure at 1 AU)
CR       = 1.8;
ASRP     = 4e-6;       % km²
muM      = 4902.8;     % km³/s²  Moon

t_span   = 30 * 86400; % 30 days, s

%% -----------------------------------------------------------------------
%% Build suite struct array
%% -----------------------------------------------------------------------
suite = struct('name',{}, 'r0',{}, 'v0',{}, 'T_orbital',{}, ...
               'e',{}, 'h_perigee',{}, 'inc',{}, 'sweep',{});
k = 0;

%% --- Group 1: Circular anchors (e=0.001, circular by convention) --------
anchors = [
    400,   0.001, 28.5;   % LEO-1
    800,   0.001, 98.0;   % LEO-2
    20200, 0.001, 55.0;   % MEO-1
    35786, 0.001,  0.1];  % GEO-1
anchor_names = {'LEO-1','LEO-2','MEO-1','GEO-1'};

for j = 1:4
    k    = k + 1;
    h    = anchors(j,1);
    e    = anchors(j,2);
    i    = anchors(j,3);
    a    = Re + h;
    [r0, v0] = posnvelos(a, e, i*pi/180, 0, 0, 0, mu);
    suite(k).name      = anchor_names{j};
    suite(k).r0        = r0;
    suite(k).v0        = v0;
    suite(k).T_orbital = 2*pi * sqrt(a^3 / mu);
    suite(k).e         = e;
    suite(k).h_perigee = h;
    suite(k).inc       = i;
    suite(k).sweep     = 'anchor';
end

%% --- Group 2: HEO eccentricity sweep (perigee=500 km, i=28.5°) ---------
% Apogee altitudes (km); HEO-E10 breaks the sweep (Molniya reference)
heo_e_params = [
    500,  2000, 28.5;   % E01
    500,  5000, 28.5;   % E02
    500, 10000, 28.5;   % E03
    500, 15000, 28.5;   % E04
    500, 20000, 28.5;   % E05
    500, 25000, 28.5;   % E06
    500, 35786, 28.5;   % E07  (anchor between E-sweep and P-sweep)
    500, 50000, 28.5;   % E08
    500, 70000, 28.5;   % E09
    300, 39000, 63.4];  % E10  Molniya reference

for j = 1:10
    k    = k + 1;
    hp   = heo_e_params(j,1);
    ha   = heo_e_params(j,2);
    i    = heo_e_params(j,3);
    rp   = Re + hp;
    ra   = Re + ha;
    a    = (rp + ra) / 2;
    e    = (ra - rp) / (ra + rp);
    [r0, v0] = posnvelos(a, e, i*pi/180, 0, 0, 0, mu);
    suite(k).name      = sprintf('HEO-E%02d', j);
    suite(k).r0        = r0;
    suite(k).v0        = v0;
    suite(k).T_orbital = 2*pi * sqrt(a^3 / mu);
    suite(k).e         = e;
    suite(k).h_perigee = hp;
    suite(k).inc       = i;
    suite(k).sweep     = 'ecc';
end
suite(k).sweep = 'molniya';  % E10 is Molniya reference, mark separately

%% --- Group 3: HEO perigee sweep (apogee=35786 km, i=28.5°) -------------
heo_p_params = [
     200, 35786, 28.5;   % P01
     400, 35786, 28.5;   % P02  (similar to HEO-E07 perigee=500 — intentionally close)
     800, 35786, 28.5;   % P03
    1200, 35786, 28.5;   % P04
    2000, 35786, 28.5;   % P05
    4000, 35786, 28.5;   % P06
    8000, 35786, 28.5;   % P07
   15000, 35786, 28.5];  % P08

for j = 1:8
    k    = k + 1;
    hp   = heo_p_params(j,1);
    ha   = heo_p_params(j,2);
    i    = heo_p_params(j,3);
    rp   = Re + hp;
    ra   = Re + ha;
    a    = (rp + ra) / 2;
    e    = (ra - rp) / (ra + rp);
    [r0, v0] = posnvelos(a, e, i*pi/180, 0, 0, 0, mu);
    suite(k).name      = sprintf('HEO-P%02d', j);
    suite(k).r0        = r0;
    suite(k).v0        = v0;
    suite(k).T_orbital = 2*pi * sqrt(a^3 / mu);
    suite(k).e         = e;
    suite(k).h_perigee = hp;
    suite(k).inc       = i;
    suite(k).sweep     = 'perigee';
end

%% --- Group 4: HEO inclination sweep (perigee=500 km, apogee=35786 km) --
% HEO-I03 uses i=35° (not 28.5°) to make all 30 orbits geometrically unique.
% (HEO-E07 already covers perigee=500/apogee=35786/i=28.5°.)
heo_i_params = [
     0.0;   % I01
     7.0;   % I02
    35.0;   % I03  (28.5° would duplicate HEO-E07 — changed to 35°)
    45.0;   % I04
    63.4;   % I05
    75.0;   % I06
    90.0;   % I07
    98.0];  % I08

for j = 1:8
    k    = k + 1;
    hp   = 500;
    ha   = 35786;
    i    = heo_i_params(j);
    rp   = Re + hp;
    ra   = Re + ha;
    a    = (rp + ra) / 2;
    e    = (ra - rp) / (ra + rp);
    [r0, v0] = posnvelos(a, e, i*pi/180, 0, 0, 0, mu);
    suite(k).name      = sprintf('HEO-I%02d', j);
    suite(k).r0        = r0;
    suite(k).v0        = v0;
    suite(k).T_orbital = 2*pi * sqrt(a^3 / mu);
    suite(k).e         = e;
    suite(k).h_perigee = hp;
    suite(k).inc       = i;
    suite(k).sweep     = 'inc';
end

assert(k == 30, 'Expected 30 suite elements, got %d', k);

save('benchmark_suite.mat', 'suite');
fprintf('Saved benchmark_suite.mat (%d orbits)\n', numel(suite));

%% -----------------------------------------------------------------------
%% Generate reference trajectories
%% -----------------------------------------------------------------------
if ~exist('ref_trajectories', 'dir')
    mkdir('ref_trajectories');
end

force_model = 'J6dragSRPmoon';

f_ref = @(t,x) orbit_eq_J6_drag_SRP_moon(t, x, mu, CD_drag, A_drag, m_drag, ...
                   Re, J, jdepoch, rhoSRP, CR, ASRP, muM);

opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

n_total = numel(suite);
summary_names  = cell(n_total, 1);
summary_tspan  = zeros(n_total, 1);
summary_rnorm  = zeros(n_total, 1);
summary_vnorm  = zeros(n_total, 1);
summary_status = cell(n_total, 1);

for k = 1:numel(suite)
    X0 = [suite(k).r0, suite(k).v0];   % 1×6
    dt = suite(k).T_orbital / 2048;
    t_eval = (0 : dt : t_span)';
    if t_eval(end) < t_span
        t_eval(end+1) = t_span;
    end

    fprintf('[%2d/%2d] %-10s ... ', k, n_total, suite(k).name);
    t_tic = tic;
    [t_ref, X_ref] = ode78(f_ref, t_eval, X0, opts);
    elapsed = toc(t_tic);
    fprintf('%.1f s\n', elapsed);

    orbit = suite(k);
    fname = fullfile('ref_trajectories', sprintf('ref_%s.mat', suite(k).name));
    save(fname, 't_ref', 'X_ref', 'orbit', 'force_model');

    summary_names{k}  = suite(k).name;
    summary_tspan(k)  = t_span / 86400;
    summary_rnorm(k)  = norm(X_ref(end, 1:3));
    summary_vnorm(k)  = norm(X_ref(end, 4:6));
    failed = any(isnan(X_ref(end,:)) | isinf(X_ref(end,:)));
    summary_status{k} = 'OK';
    if failed
        summary_status{k} = 'FAILED';
    end
end

%% Summary table
fprintf('\n%-10s | t_span(d) | |r_final|(km) | |v_final|(km/s) | Status\n', 'Name');
fprintf('%s\n', repmat('-', 1, 62));
n_failed = 0;
for k = 1:n_total
    fprintf('%-10s | %9.1f | %12.3f | %15.6f | %s\n', ...
        summary_names{k}, summary_tspan(k), summary_rnorm(k), ...
        summary_vnorm(k), summary_status{k});
    if strcmp(summary_status{k}, 'FAILED')
        n_failed = n_failed + 1;
    end
end
fprintf('%s\n', repmat('-', 1, 62));
if n_failed == 0
    fprintf('All %d reference trajectories OK.\n', n_total);
else
    fprintf('WARNING: %d / %d trajectories FAILED.\n', n_failed, n_total);
end
