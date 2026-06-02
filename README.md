# MGPI — Modified Gauss-Picard Integration

> Companion code for:  
> E. Denenberg, **"Legendre Basis Functions Improve Picard Iteration Accuracy for Highly Eccentric Orbit Propagation"**  
> *AIAA Journal of Guidance, Control, and Dynamics* (submitted June 2026)

---

## Overview

MGPI generalises the Modified Picard-Chebyshev Integration (MCPI) method to arbitrary Gegenbauer polynomial bases via a single blending parameter **α**.

- **α = 0** recovers standard MCPI with Chebyshev-Gauss-Lobatto (CGL) collocation nodes  
- **α > 0** shifts nodes toward the interior of each segment (Gegenbauer-Gauss-Lobatto, GGL), improving convergence for highly eccentric orbits

This repository contains the full benchmark suite used in the paper:

- **30 orbital scenarios** — 4 circular anchors (LEO×2, MEO, GEO) + 10-orbit eccentricity sweep + 8-orbit perigee sweep + 8-orbit inclination sweep
- **Force model** — J2–J6 zonal harmonics + atmospheric drag + solar radiation pressure + lunar third body (`orbit_eq_J6_drag_SRP_moon.m`)
- **Propagation span** — 30 days per orbit (alpha and sensitivity sweeps), 7 days (full timing sweep)
- **Reference accuracy** — ode78 at RelTol = AbsTol = 1×10⁻¹³

---

## Requirements

- **MATLAB R2020b or later**  
  The Curve Fitting Toolbox is required for `ode78`, which is used to generate reference trajectories.
- **Python 3.9 or later** (standard library only)  
  Required only for the LaTeX figure post-processing step (`gen_figures_tex.py`).

---

## Dependencies

The following libraries must be downloaded and placed in accessible folders.  
Edit `setup_matlab_paths.m` to point MATLAB to each one.

| Library | Where to get it | What it is used for |
|---------|----------------|---------------------|
| **Vallado's Astrodynamics Library** | Clone `https://github.com/CelesTrak/fundamentals-of-astrodynamics` and sparse-checkout `software/matlab/` | `getgravc` (WGS-84 constants), `ode78` (reference integrator) |
| **OPQ** (Gautschi) | Download from https://people.math.ethz.ch/~waldvoge/Programs/Gautschi/ | Orthogonal polynomial recurrence + Gauss-Lobatto quadrature nodes for GGL collocation |
| **matlab2tikz** | Clone `https://github.com/matlab2tikz/matlab2tikz` | Export MATLAB figures to PGFPlots/TikZ |
| **table2latex** | Search MATLAB File Exchange for *"table2latex"* | Export MATLAB tables to LaTeX `tabular` |
| **maxdistcolor** | Search MATLAB File Exchange for *"maxdistcolor"* (S. Cobeldick) | Perceptually maximally distinct color palette for α-sweep figures |

### Sparse-checkout for Vallado's library (saves bandwidth)

```bash
git clone --no-checkout https://github.com/CelesTrak/fundamentals-of-astrodynamics.git
cd fundamentals-of-astrodynamics
git sparse-checkout init
git sparse-checkout set software/matlab
git checkout
```

---

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/BraudeMechInd4-0/MGPI.git
   cd MGPI
   ```

2. Download/clone all five dependencies listed above.

3. Open `setup_matlab_paths.m` and update the five path variables to your local copies.

4. In MATLAB:
   ```matlab
   run setup_matlab_paths
   ```

---

## Quick start (5–10 min)

```matlab
run setup_matlab_paths
sanity_check        % verifies integrator output structure on a simple case
test_sweep_alpha    % runs 2 cases (LEO-1, α = {0, 0.5}) as a smoke test
```

---

## Running the benchmark

### Full pipeline (one command, ~1–4 days)

```matlab
run_all_sweeps
```

### Stage by stage

| Stage | Script | Runtime | Output |
|-------|--------|---------|--------|
| 0 — Reference trajectories | `init_benchmark_suite` | ~1–2 h | `ref_trajectories/` + `benchmark_suite.mat` |
| 1 — Alpha sweep | `sweep_alpha` | ~1–2 h | `results/sweep/` |
| 2 — Sensitivity sweep | `sweep_sensitivity` | ~1–2 h | `results/sensitivity/` |
| 3 — Full timing sweep | `sweep_full` | ~6–12 h | `results/full/` |
| 4 — Figures & tables | `regen_figures` | ~30 min | `results/figures/` |

---

## Generating paper figures

After Stage 4:

```bash
python gen_figures_tex.py        # assembles results/figures_section.tex
```

Copy `results/figures_processed/` to the `Figures/` folder of your LaTeX project and `\input` `figures_section.tex`.

For subplot panels (e.g., 2×2 full-sweep comparison):

```bash
python make_subplot_tikz.py E03 E05 E07 E09
```

---

## `odeMPGI` API

```matlab
[tout, xout, err, model] = odeMPGI(ODEFUN, TSPAN, X0)
[tout, xout, err, model] = odeMPGI(ODEFUN, TSPAN, X0, options)
```

**Inputs**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `ODEFUN` | Function handle `f(t,x)` for `dx/dt = f` | — |
| `TSPAN` | `[t0, tf]` or vector of output times | — |
| `X0` | Initial state vector | — |
| `options.alpha` | Gegenbauer parameter (0 = MCPI/Chebyshev) | `0.5` |
| `options.N` | Nodes per segment | `16` |
| `options.deltaT` | Segment duration (s) | `T_orbital / 16` |
| `options.AbsTol` | Picard convergence absolute tolerance | `1e-12` |
| `options.RelTol` | Picard convergence relative tolerance | `1e-9` |

**Outputs**

| Variable | Description |
|----------|-------------|
| `tout` | Output time vector (s) |
| `xout` | State history — rows correspond to `tout` |
| `err` | Error code: `0` = converged, `-1` = failed |
| `model` | Struct with segment timestamps (`model.ts`) and polynomial coefficients (`model.b`) |

**Minimal example — LEO propagation**

```matlab
run setup_matlab_paths

mu   = 398600.4418;          % km^3/s^2
r0   = [6778; 0; 0];         % km, 400 km altitude
v0   = [0; sqrt(mu/norm(r0)); 0];  % km/s, circular
X0   = [r0; v0];
T    = 2*pi*sqrt(norm(r0)^3/mu);   % orbital period (s)

opts.alpha  = 0.5;
opts.N      = 16;
opts.deltaT = T / 4;         % 4 segments per orbit

[t, x] = odeMPGI(@orbit_eq, [0, 10*T], X0, opts);

plot(x(:,1), x(:,2))
xlabel('x (km)'); ylabel('y (km)'); axis equal
```

---

## Repository structure

```
MGPI_forJGCD/
├── README.md
├── setup_matlab_paths.m
├── run_all_sweeps.m          pipeline orchestrator
├── regen_figures.m           Stage 4: regenerate all figures
├── init_benchmark_suite.m    Stage 0: reference trajectories
├── odeMPGI.m                 main integrator
├── odeMPCI.m                 MCPI (α = 0 special case)
├── orbit_eq*.m               force models
├── drag_accel.m  moon_v.m  sun_v.m
├── kep_elements.m  posnvelos.m
├── gegenbauer_precompute.m  gegenbauer_eval.m
├── sweep_alpha.m  sweep_sensitivity.m  sweep_full.m
├── run_sweep_case.m
├── plot_alpha_sweep.m  plot_sensitivity.m  plot_full.m
├── get_alpha_colors.m
├── sanity_check.m  test_sweep_alpha.m
├── gen_figures_tex.py
└── make_subplot_tikz.py
```

---

## License

TBD