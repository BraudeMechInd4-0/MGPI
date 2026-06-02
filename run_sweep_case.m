function [gbl, seg, nod] = run_sweep_case(orb, fm, alpha_val, delta, N_nodes, t_ref, X_ref, C)
% One (orbit, force_model, alpha, delta, N) run.
% orb     : suite element (fields: r0, v0, T_orbital)
% fm      : 'J6' or 'J6drag'
% alpha_val : Gegenbauer alpha (0 = Chebyshev)
% delta   : number of segments per orbit
% N_nodes : nodes per segment (options.N)
% t_ref   : reference time vector (L×1 s)
% X_ref   : reference trajectory (L×6 km and km/s)
% C       : constants struct (mu, Re, J, CD_drag, A_drag, m_drag)
%           C.debug = true  →  store seg.f_osc, seg.evec and plot orbital
%                              elements vs time for inspection

if strcmp(fm, 'J6')
    f_ode = @(t,x) orbit_eq_J6_drag(t, x, C.mu, 0,         C.A_drag, C.m_drag, C.Re, C.J);
elseif strcmp(fm, 'Kep')
    f_ode = @(t,x) orbit_eq(t, x, C.mu);
elseif strcmp(fm, 'J6dragSRPmoon')
    f_ode = @(t,x) orbit_eq_J6_drag_SRP_moon(t, x, C.mu, C.CD_drag, C.A_drag, C.m_drag, ...
                       C.Re, C.J, C.jdepoch, C.rhoSRP, C.CR, C.ASRP, C.muM);
else
    f_ode = @(t,x) orbit_eq_J6_drag(t, x, C.mu, C.CD_drag, C.A_drag, C.m_drag, C.Re, C.J);
end

opts.alpha   = alpha_val;
opts.N       = N_nodes;
opts.Sec     = orb.T_orbital / delta;
opts.AbsTol  = 1e-9;
opts.RelTol  = 1e-6;
opts.maxIter = 2000;

t_tic = tic;
try
    [tout, xout, errCode, model] = odeMPGI(f_ode, t_ref, [orb.r0, orb.v0], opts);
catch
    [gbl, seg, nod] = failed_result(alpha_val, N_nodes, -3);
    return
end
wt = toc(t_tic);

% Trim NaN tail
valid = ~isnan(tout);
tout  = tout(valid);
xout  = xout(valid, :);
M     = numel(model.iters);

% Satellite decay: r < 0.5*Re means the integrator wandered inside Earth
if any(sqrt(sum(xout(:,1:3).^2, 2)) < 0.5*C.Re)
    [gbl, seg, nod] = failed_result(alpha_val, N_nodes, -3);
    return
end

%% Global metrics (xout and X_ref both at t_ref times — direct subtraction)
err_r = sqrt(sum((xout(:,1:3) - X_ref(1:numel(tout),1:3)).^2, 2));   % L×1
err_v = sqrt(sum((xout(:,4:6) - X_ref(1:numel(tout),4:6)).^2, 2));   % L×1

gbl.alpha        = alpha_val;
gbl.r_err_final  = err_r(end);
gbl.v_err_final  = err_v(end);
gbl.total_iters  = sum(model.iters);
gbl.total_fevals = gbl.total_iters * N_nodes;
gbl.wall_time    = wt;
gbl.error_code   = errCode;

%% Per-segment metrics
seg_idx = max(1, min(M, floor((tout - model.ts(1)) / opts.Sec) + 1));  % L×1

seg.r_err  = accumarray(seg_idx, err_r, [M,1], @max);
seg.v_err  = accumarray(seg_idx, err_v, [M,1], @max);
seg.iters  = model.iters;
seg.fevals = model.iters * N_nodes;

% Orbital state at segment start from actual solution (first point of each segment)
t_seg_start = model.ts(1:M);                                        % M×1
seg_first_cell = arrayfun(@(k) find(seg_idx == k, 1, 'first'), (1:M)', 'UniformOutput', false);
if any(cellfun(@isempty, seg_first_cell))
    [gbl, seg, nod] = failed_result(alpha_val, N_nodes, -3);
    return
end
seg_first = cell2mat(seg_first_cell);
rv = xout(seg_first, 1:3);
vv = xout(seg_first, 4:6);
r_norm = sqrt(sum(rv.^2, 2));    % ‖r‖  M×1
v_norm = sqrt(sum(vv.^2, 2));    % ‖v‖  M×1

% Vectorised orbital elements — mirrors kep_elements.m lines 21, 32, 52-54
hv    = cross(rv, vv, 2);                                           % h = cross(r,v)  [kep_elements line 21]
evec  = cross(vv, hv, 2) ./ C.mu - rv ./ r_norm;                   % ev  [kep_elements line 32]
e_osc = sqrt(sum(evec.^2, 2));

rdotv = sum(rv .* vv, 2);

% E3: Kepler's equation  M = E − e·sin(E), vectorised Newton from M = n·t
n_mean = 2*pi / orb.T_orbital;
M_osc  = mod(n_mean .* t_seg_start, 2*pi);
E_iter = M_osc + e_osc .* sin(M_osc);
for iter = 1:10
    dE     = (M_osc - E_iter + e_osc .* sin(E_iter)) ./ (1 - e_osc .* cos(E_iter));
    E_iter = E_iter + dE;
    if max(abs(dE)) < 1e-14, break; end
end
E3_osc = mod(E_iter + pi, 2*pi) - pi;

E_osc = E3_osc;

% f_osc for debug plot (eccentricity vector direction)
cosf  = sum(evec .* rv, 2) ./ max(e_osc .* r_norm, 1e-12);
cosf  = max(-1, min(1, cosf));
f_osc = acos(cosf);
f_osc(rdotv < 0) = 2*pi - f_osc(rdotv < 0);

seg.t_start  = t_seg_start;
seg.r_start  = r_norm;
seg.v_start  = v_norm;
seg.rv_start = rv;          % M×3 position vectors (km)
seg.vv_start = vv;          % M×3 velocity vectors (km/s)
seg.e_osc    = e_osc;
seg.f_osc    = f_osc;
seg.evec     = evec;
seg.E_osc    = E_osc;
seg.E3_osc   = E3_osc;

debug_E = isfield(C, 'debug_E') && C.debug_E;
if debug_E
    % E1: orbit-equation formula  cos(f) = (p/r - 1)/e  where p = h²/μ
    p_osc  = sum(hv.^2, 2) ./ C.mu;
    cosf1  = max(-1, min(1, (p_osc ./ r_norm - 1) ./ max(e_osc, 1e-12)));
    f1_osc = acos(cosf1);
    f1_osc(rdotv < 0) = 2*pi - f1_osc(rdotv < 0);
    seg.E1_osc = mod(atan2(sqrt(max(0,1-e_osc.^2)).*sin(f1_osc), e_osc+cos(f1_osc)) + pi, 2*pi) - pi;

    % E2: direct vis-viva  E = atan2(r·v·sqrt(a/μ), a−r)
    a_sma      = 1 ./ (2./r_norm - v_norm.^2 ./ C.mu);
    seg.E2_osc = atan2(rdotv .* sqrt(a_sma ./ C.mu), a_sma - r_norm);
end

%% Per-node metrics (xout already evaluated at t_ref points — no reconstruction needed)
nod.r_err   = err_r;
nod.v_err   = err_v;
nod.t       = tout;
nod.seg_idx = seg_idx;

%% Debug plots (only when C.debug == true)
if isfield(C, 'debug') && C.debug
    t_days = t_seg_start / 86400;
    ttl = sprintf('%s / %s / \alpha=%.2f', orb.name, fm, alpha_val);

    figure('Name', ['Elements: ' ttl], 'Position', [50 50 1400 900]);

    subplot(3,2,1);
    plot(t_days, e_osc, '.-');
    xlabel('t (days)'); ylabel('e_{osc}');
    title('Osculating eccentricity'); grid on;

    subplot(3,2,2);
    plot(t_days, rad2deg(f_osc), '.-');
    xlabel('t (days)'); ylabel('f_{osc} (deg)');
    title('True anomaly'); grid on;

    subplot(3,2,3);
    plot(t_days, rad2deg(E_osc), '.-');
    xlabel('t (days)'); ylabel('E_{osc} (deg)');
    title('Eccentric anomaly (raw, before unwrap)'); grid on;

    subplot(3,2,4);
    plot(t_days, rad2deg(unwrap(E_osc)), '.-');
    xlabel('t (days)'); ylabel('E_{osc} unwrapped (deg)');
    title('Eccentric anomaly (unwrapped)'); grid on;

    subplot(3,2,5);
    plot(t_days, evec(:,1), '.-', t_days, evec(:,2), '.-', t_days, evec(:,3), '.-');
    xlabel('t (days)'); ylabel('evec component');
    legend('e_x','e_y','e_z'); title('Eccentricity vector components'); grid on;

    subplot(3,2,6);
    plot(t_days, rdotv, '.-');
    xlabel('t (days)'); ylabel('r \cdot v (km^2/s)');
    title('r\cdotv (sign determines f quadrant)'); grid on;

    sgtitle(ttl);
end

end % run_sweep_case

% -------------------------------------------------------------------------
function [gbl, seg, nod] = failed_result(alpha_val, ~, code)
gbl.alpha        = alpha_val;
gbl.r_err_final  = NaN;
gbl.v_err_final  = NaN;
gbl.total_iters  = NaN;
gbl.total_fevals = NaN;
gbl.wall_time    = NaN;
gbl.error_code   = code;
seg = struct();
nod = struct();
end
