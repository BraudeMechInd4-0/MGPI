#!/usr/bin/env python3
"""
gen_figures_tex.py

Generates results/figures_section.tex - a LaTeX fragment (\\input-able, no preamble)
with three \\subsection blocks of figures from the MGPI sweep results.

Also writes preprocessed tikz files to results/figures_processed/:
  - title commented out (text reused as caption)
  - \addlegendentry lines commented out (legend moved to a single legend figure per section)

Usage:
    python gen_figures_tex.py
Then copy results/figures_processed/ contents to Figures/ in your LaTeX project.
"""

import re
import shutil
from pathlib import Path

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
PROCESS_TIKZ = False        # True: preprocess tikz files and copy to DST_ROOT
                            # False: only generate the .tex file
TABLE_SIGFIGS = 3           # Significant figures for \num{} values in tables (None = no rounding)

MARKER_SCALE = {            # Per-shape multipliers (sensitivity + full tikz files)
    '*':         2.5,       # filled circle — baseline
    'square*':   2.5,       # filled square
    'triangle*': 4.0,       # filled triangle — boost more (appears small)
    'diamond*':  1.5,       # filled diamond — boost less (appears large)
    'square':    2.5,       # open square  (ode45 ref in full)
    'triangle':  4.0,       # open triangle (ode45 ref in full)
    'diamond':   1.5,       # open diamond  (ode45 ref in full)
    'x':         2.5,       # cross (MPCM failed in full)
    'default':   2.5,
}
SENS_MARKER_SCALE = MARKER_SCALE   # alias kept for clarity

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SRC_ROOT          = Path('results/figures')
DST_ROOT          = Path('results/figures_processed')
TEX_OUT           = Path('results/figures_section.tex')

SRC_ERROR_TIME    = SRC_ROOT / 'alpha_sweep' / 'error_time'
SRC_RUNTIME       = SRC_ROOT / 'alpha_sweep' / 'runtime'
SRC_LEGEND_ALPHA  = SRC_ROOT / 'alpha_sweep' / 'legend.tikz'
SRC_SENSITIVITY   = SRC_ROOT / 'sensitivity'
SRC_FULL          = SRC_ROOT / 'full'

DST_ERROR_TIME    = DST_ROOT / 'error_time'
DST_RUNTIME       = DST_ROOT / 'runtime'
DST_SENSITIVITY   = DST_ROOT / 'sensitivity'
DST_FULL          = DST_ROOT / 'full'
SRC_TABLES         = Path('results/tables')
DST_TABLES         = DST_ROOT / 'Tables'

DST_LEGEND_ALPHA   = DST_ROOT / 'legend_alpha.tikz'
DST_LEGEND_RUNTIME = DST_ROOT / 'legend_runtime.tikz'
DST_LEGEND_SENS    = DST_ROOT / 'legend_sensitivity.tikz'
DST_LEGEND_FULL    = DST_ROOT / 'legend_full.tikz'

# Orbits to include in figures_section.tex for the full-sweep section.
# All full/*.tikz files are always processed; only these appear in the tex.
FULL_TEX_ORBITS = ['HEO-E07', 'HEO-E08']

# ---------------------------------------------------------------------------
# Tikz preprocessing
# ---------------------------------------------------------------------------

def extract_title(content: str) -> str:
    m = re.search(r'\btitle=\{', content)
    if not m:
        return ''
    start = m.end()
    depth = 1
    i = start
    while i < len(content) and depth > 0:
        if content[i] == '{':
            depth += 1
        elif content[i] == '}':
            depth -= 1
        i += 1
    return content[start:i-1].strip() if depth == 0 else ''


def round_num_macros(content: str, sigfigs: int) -> str:
    """Round every \\num{X} value in a LaTeX table to sigfigs significant figures."""
    def _replace(m):
        try:
            val = float(m.group(1))
            return r'\num{' + f'{val:.{sigfigs - 1}e}' + '}'
        except ValueError:
            return m.group(0)
    return re.sub(r'\\num\{([^}]+)\}', _replace, content)


def copy_table(src: Path, dst: Path):
    """Copy a table tex file to dst, applying TABLE_SIGFIGS rounding if set."""
    content = src.read_text(encoding='utf-8')
    if TABLE_SIGFIGS is not None:
        content = round_num_macros(content, TABLE_SIGFIGS)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(content, encoding='utf-8')


def scale_marker_sizes(content: str, scale) -> str:
    """Multiply every mark size=X.XXXXpt value by scale.
    scale may be a float (uniform) or a dict {mark_name: factor, 'default': factor}.
    """
    if isinstance(scale, (int, float)):
        def _replace(m):
            return f'mark size={float(m.group(1)) * scale:.4f}pt'
        return re.sub(r'mark size=([\d.]+)pt', _replace, content)

    # Per-shape: inspect each line for mark= to pick the right factor
    lines = content.splitlines(keepends=True)
    out = []
    for line in lines:
        if 'mark size=' in line and r'\addplot' in line:
            mark_m = re.search(r'\bmark=([\w*]+)', line)
            mark_name = mark_m.group(1) if mark_m else 'default'
            s = scale.get(mark_name, scale.get('default', 1.0))
            def _sz(m, s=s):
                return f'mark size={float(m.group(1)) * s:.4f}pt'
            line = re.sub(r'mark size=([\d.]+)pt', _sz, line)
        out.append(line)
    return ''.join(out)



def patch_anchor_tikz(src: Path, dst: Path) -> str:
    """
    Fix anchor runtime tikz files: spread orbit groups wider and tighten per-alpha offsets.

    Old layout: group centers at 1,2,3,4  offset = (j-6)*0.14  (j=1..11)
    New layout: group centers at 1,4,7,10 offset = (j-6)*0.07

    Transform: for each x value, recover group rank (0-3) and offset, then apply new formula.
    Also updates xmin/xmax/xtick axis parameters.
    """
    content = src.read_text(encoding='utf-8')
    title = extract_title(content)

    # Transform x-coordinates in each addplot table block
    def transform_series(block: str) -> str:
        rows = re.findall(r'^([\d.\-]+)\t([\d.eE+\-]+)\\\\$', block, re.MULTILINE)
        if not rows:
            return block
        xs = [float(r[0]) for r in rows]
        ys = [r[1] for r in rows]
        # Sort by x to determine group rank (0=lowest x → group 1, etc.)
        order = sorted(range(len(xs)), key=lambda i: xs[i])
        offset = xs[order[0]] - 1.0          # offset from group-1 center
        new_centers = [1, 4, 7, 10]
        new_offset  = offset * (18.0 / 14.0)  # scale offset: old mult=0.14, new mult=0.18
        new_xs = [None] * len(xs)
        for rank, orig_i in enumerate(order):
            new_xs[orig_i] = new_centers[rank] + new_offset
        new_rows = '\n'.join(f'{x:.6g}\t{y}\\\\' for x, y in zip(new_xs, ys)) + '\n'
        return re.sub(r'([\d.\-]+\t[\d.eE+\-]+\\\\[\n\r])+',
                      lambda m: new_rows, block)

    # Apply to each table block
    content = re.sub(
        r'(table\[row sep=crcr\]\{%\n)(.*?)(\};)',
        lambda m: m.group(1) + transform_series(m.group(2)) + m.group(3),
        content, flags=re.DOTALL
    )

    # Update axis x range and ticks
    content = re.sub(r'xmin=[\d.\-]+,', 'xmin=-0.5,', content)
    content = re.sub(r'xmax=[\d.\-]+,', 'xmax=11.5,', content)
    content = re.sub(r'xtick=\{[^}]+\}', 'xtick={1,4,7,10}', content)

    # Run standard preprocess (title, width, definecolors)
    deferred_colors = []
    out = []
    for line in content.splitlines(keepends=True):
        if re.match(r'\\definecolor', line.rstrip('%').strip()):
            deferred_colors.append(line)
        elif re.match(r'\\begin\{tikzpicture\}', line.strip()):
            out.append(line); out.extend(deferred_colors); deferred_colors = []
        elif re.match(r'\s*title style\s*=', line):
            out.append('%' + line)
        elif re.match(r'\s*title\s*=\s*\{', line):
            out.append('%' + line)
        elif re.match(r'\s*\\addlegendentry', line):
            out.append('%' + line)
        elif re.match(r'\s*at\s*=\s*\{', line):
            out.append('%' + line)
        elif re.match(r'\s*width\s*=', line):
            indent = re.match(r'(\s*)', line).group(1)
            out.append(f'{indent}width=\\linewidth,\n')
        else:
            out.append(line)

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(''.join(out), encoding='utf-8')
    return title


def preprocess_tikz(src: Path, dst: Path, strip_legend: bool = True,
                    marker_scale: float = None,
                    width: str = None, height: str = None) -> str:
    """
    Comment out title style/value lines (always) and optionally \\addlegendentry lines.
    Returns the extracted title text for use as a caption.
    strip_legend=False keeps legend entries intact (used when no legend figure is provided).
    marker_scale: if given, all mark size= values are multiplied by this factor.
    width: override replacement width string (default: \\linewidth).
    height: if given, replace height= line with this value; otherwise keep original.
    """
    content = src.read_text(encoding='utf-8')
    title = extract_title(content)

    deferred_colors = []
    out = []
    for line in content.splitlines(keepends=True):
        if re.match(r'\\definecolor', line.rstrip('%').strip()):
            deferred_colors.append(line)
        elif re.match(r'\\begin\{tikzpicture\}', line.strip()):
            out.append(line)
            out.extend(deferred_colors)
            deferred_colors = []
        elif re.match(r'\s*title style\s*=', line):
            out.append('%' + line)
        elif re.match(r'\s*title\s*=\s*\{', line):
            out.append('%' + line)
        elif strip_legend and re.match(r'\s*\\addlegendentry', line):
            out.append('%' + line)
        elif re.match(r'\s*at\s*=\s*\{', line):
            out.append('%' + line)
        elif re.match(r'\s*width\s*=', line):
            indent = re.match(r'(\s*)', line).group(1)
            w = width if width is not None else r'\linewidth'
            out.append(f'{indent}width={w},\n')
        elif re.match(r'\s*height\s*=', line):
            if height is not None:
                indent = re.match(r'(\s*)', line).group(1)
                out.append(f'{indent}height={height},\n')
            else:
                out.append(line)
        else:
            out.append(line)

    result = ''.join(out)
    if marker_scale is not None:
        result = scale_marker_sizes(result, marker_scale)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(result, encoding='utf-8')
    return title


def build_alpha_legend(src_data_tikz: Path, src_legend_tikz: Path, dst: Path):
    """
    Build a proper legend for the alpha sweep.

    MATLAB's legend.tikz only has \\addlegendentry lines (no \\addplot), which
    breaks pgfplots at \\end{axis}.  Fix: extract plot styles (including their
    \\definecolor definitions) from one data tikz (where all plots are
    forget-plot) and pair them with the \\addlegendentry texts from legend.tikz.
    """
    data_content = src_data_tikz.read_text(encoding='utf-8')
    leg_content  = src_legend_tikz.read_text(encoding='utf-8')

    # Collect \definecolor lines (they sit outside \begin{tikzpicture})
    color_defs = [ln for ln in data_content.splitlines()
                  if re.match(r'\\definecolor', ln.rstrip('%').strip())]

    # Extract unique plot styles in order of first appearance (strip forget plot)
    seen_styles = []
    seen_set = set()
    for line in data_content.splitlines():
        if re.match(r'\s*\\addplot', line) and 'forget plot' in line:
            style = line.strip().rstrip(';')
            style = re.sub(r',?\s*forget plot', '', style)
            style = re.sub(r'\s*coordinates\s*\{.*$', '', style).strip()
            if style not in seen_set:
                seen_set.add(style)
                seen_styles.append(style)

    # Extract \addlegendentry texts from legend.tikz
    leg_texts = []
    for line in leg_content.splitlines():
        m = re.search(r'\\addlegendentry\{(.*)\}', line)
        if m:
            leg_texts.append(m.group(1))

    pairs = list(zip(seen_styles, leg_texts))
    if not pairs:
        print(f'  WARNING: build_alpha_legend found no pairs.')

    out = [
        '% Auto-generated alpha legend (gen_figures_tex.py)',
        r'\begin{tikzpicture}',
    ]
    out += color_defs
    out += [
        r'\begin{axis}[',
        r'  hide axis,',
        r'  scale only axis,',
        r'  width=0.001pt, height=0.001pt,',
        r'  xmin=0, xmax=1, ymin=0, ymax=1,',
        r'  legend columns=-1,',
        r'  legend style={at={(0,0)}, anchor=south west, draw=none, /tikz/every even column/.append style={column sep=0.4cm}}',
        r']',
    ]
    for style, leg_text in pairs:
        out.append(f'  {style} coordinates {{(0,0)}};')
        out.append(f'  \\addlegendentry{{{leg_text}}}')
    out += [r'\end{axis}', r'\end{tikzpicture}', '']

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text('\n'.join(out), encoding='utf-8')
    print(f'  Built alpha legend ({len(pairs)} entries) -> {dst}')


def build_legend(src_tikz: Path, dst: Path, label: str = 'legend',
                 marker_scale: float = None):
    """
    Generate a legend-only tikzpicture by scanning src_tikz for
    \\addplot [options] / \\addlegendentry pairs.

    Tracks the most recent non-forget-plot \\addplot style line as we scan
    forward, so coordinate blocks of any length are handled correctly.
    Only the options line (first line of the addplot block) is kept;
    coordinates are replaced with {}.
    """
    content = src_tikz.read_text(encoding='utf-8')
    lines = content.splitlines()

    color_defs = [ln for ln in lines
                  if re.match(r'\\definecolor', ln.rstrip('%').strip())]

    pairs = []
    last_plot_style = None
    for line in lines:
        if re.match(r'\s*\\addplot', line) and 'forget plot' not in line:
            # Keep only the options — strip data format (coordinates/table) and any trailing {
            style = line.strip().rstrip(';')
            style = re.sub(r'\s*(?:coordinates|table)\b.*$', '', style)
            last_plot_style = style.strip()
        elif re.match(r'\s*\\addlegendentry', line):
            leg_m = re.search(r'\\addlegendentry\{(.*)\}', line)
            if leg_m and last_plot_style:
                pairs.append((last_plot_style + ' coordinates {(0,0)};', leg_m.group(1)))
                last_plot_style = None  # consume

    if not pairs:
        print(f'  WARNING: no legend pairs found in {src_tikz}, writing empty legend.')

    out = [
        f'% Auto-generated {label} (gen_figures_tex.py)',
        r'\begin{tikzpicture}',
    ]
    out += color_defs
    out += [
        r'\begin{axis}[',
        r'  hide axis,',
        r'  scale only axis,',
        r'  width=0.001pt, height=0.001pt,',
        r'  xmin=0, xmax=1, ymin=0, ymax=1,',
        r'  legend columns=-1,',
        r'  legend style={at={(0,0)}, anchor=south west, draw=none, /tikz/every even column/.append style={column sep=0.4cm}}',
        r']',
    ]
    for plot_line, leg_text in pairs:
        line = f'  {plot_line}'
        if marker_scale is not None:
            line = scale_marker_sizes(line, marker_scale)
        out.append(line)
        out.append(f'  \\addlegendentry{{{leg_text}}}')
    out += [r'\end{axis}', r'\end{tikzpicture}', '']

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text('\n'.join(out), encoding='utf-8')
    print(f'  Built {label} ({len(pairs)} entries) -> {dst}')


# ---------------------------------------------------------------------------
# LaTeX helpers
# ---------------------------------------------------------------------------

def _tex_path(*parts) -> str:
    return '/'.join(str(p) for p in parts)


def write_legend_figure(f, figures_prefix: str, legend_name: str):
    path = _tex_path(figures_prefix, legend_name)
    f.write('\\begin{figure}[htbp]\n')
    f.write('  \\centering\n')
    f.write(f'  \\begin{{CacheMe}}{{tikz, add dependencies={{{path}}}}}\n')
    f.write(f'  \\input{{../{path}}}\n')
    f.write(f'  \\end{{CacheMe}}\n')
    f.write('\\end{figure}\n\n')


def write_figure(f, figures_prefix: str, tikz_name: str, caption: str, label: str):
    path = _tex_path(figures_prefix, tikz_name)
    f.write('\\begin{figure}[htbp]\n')
    f.write('  \\centering\n')
    f.write(f'  \\begin{{CacheMe}}{{tikz, add dependencies={{{path}}}}}\n')
    f.write(f'  \\input{{../{path}}}\n')
    f.write(f'  \\end{{CacheMe}}\n')
    f.write(f'  \\caption{{{caption}}}\n')
    f.write(f'  \\label{{{label}}}\n')
    f.write('\\end{figure}\n\n')


def write_table(f, fname: str, caption: str, label: str):
    f.write('\\begin{table}[htbp]\n')
    f.write('  \\centering\n')
    f.write(f'  \\caption{{{caption}}}\n')
    f.write(f'  \\label{{{label}}}\n')
    f.write(f'  \\input{{Tables/{fname}}}\n')
    f.write('\\end{table}\n\n')


def slug(name: str) -> str:
    return name.replace('-', '_')


def clean_caption(title: str) -> str:
    """Strip force-model suffix and escape LaTeX specials.
    Returns '' if stripping leaves unbalanced LaTeX (caller uses fallback)."""
    title = re.sub(r'\s*[/—]\s*J6dragSRPmoon.*$', '', title).strip()
    if title.count('$') % 2 != 0 or title.count('{') != title.count('}'):
        return ''
    title = title.replace('_', r'\_')
    if title and not title.endswith('.'):
        title += '.'
    return title


# ---------------------------------------------------------------------------
# Section 1 — position error vs time
# ---------------------------------------------------------------------------

FAMILIES = [
    ('Circular Orbits',            ['LEO-1', 'LEO-2', 'MEO-1', 'GEO-1']),
    ('HEO --- Eccentricity Sweep', [f'HEO-E{i:02d}' for i in range(1, 11)]),
    ('HEO --- Perigee Sweep',      [f'HEO-P{i:02d}' for i in range(1, 9)]),
    ('HEO --- Inclination Sweep',  [f'HEO-I{i:02d}' for i in range(1, 9)]),
]


def write_section_1(f):
    f.write('% ================================================================\n')
    f.write('\\subsection{Position Error vs.\\ Time}\n\n')
    write_legend_figure(f, 'Figures', 'legend_alpha.tikz')

    for subhead, orbits in FAMILIES:
        f.write(f'\\subsubsection{{{subhead}}}\n\n')
        for orbit in orbits:
            src = SRC_ERROR_TIME / f'{orbit}_J6dragSRPmoon.tikz'
            if not src.exists():
                print(f'  SKIP (missing): {src.name}')
                continue
            title = preprocess_tikz(src, DST_ERROR_TIME / src.name) if PROCESS_TIKZ else extract_title(src.read_text(encoding='utf-8'))
            caption = clean_caption(title) or f'Position error vs.\\ time --- {orbit}.'
            write_figure(f, 'Figures/error_time', src.name,
                         caption, f'fig:errt_{slug(orbit)}')
            print(f'  error_time/{src.name}')


# ---------------------------------------------------------------------------
# Section 2 — runtime summary
# ---------------------------------------------------------------------------

RUNTIME_GROUPS = [
    ('Wall Time', [
        'wtime_ecc', 'wtime_perigee', 'wtime_inc', 'wtime_anchors',
    ]),
    ('Function Evaluations', [
        'fevals_ecc', 'fevals_perigee', 'fevals_inc', 'fevals_anchors',
    ]),
    (r'Position Error at $t_\mathrm{end}$', [
        'rerr_ecc', 'rerr_perigee', 'rerr_inc', 'rerr_anchors',
    ]),
]

PARAM_LABEL = {
    'ecc':     'eccentricity',
    'perigee': 'perigee altitude',
    'inc':     'inclination',
    'anchors': 'circular reference orbits',
}
METRIC_LABEL = {
    'wtime':  'Wall time',
    'fevals': 'Function evaluations',
    'rerr':   r'Position error at $t_\mathrm{end}$',
}


ALPHA_TABLES = [
    ('alpha_r_err_J6dragSRPmoon.tex',  r'Position error at $t_\mathrm{end}$ (km)'),
    ('alpha_wtime_J6dragSRPmoon.tex',  'Wall time (s)'),
    ('alpha_fevals_J6dragSRPmoon.tex', 'Function evaluations'),
]


def write_section_2(f):
    f.write('% ================================================================\n')
    f.write('\\subsection{Alpha Sweep: Runtime Summary}\n\n')

    # Runtime legend
    for stems in [s for _, s in RUNTIME_GROUPS]:
        for stem in stems:
            src_leg = SRC_RUNTIME / f'{stem}_J6dragSRPmoon.tikz'
            if src_leg.exists():
                if PROCESS_TIKZ:
                    build_legend(src_leg, DST_LEGEND_RUNTIME)
                write_legend_figure(f, 'Figures', 'legend_runtime.tikz')
                break
        else:
            continue
        break
    else:
        print('  WARNING: no runtime tikz found for legend.')

    for subhead, stems in RUNTIME_GROUPS:
        f.write(f'\\subsubsection{{{subhead}}}\n\n')
        for stem in stems:
            src = SRC_RUNTIME / f'{stem}_J6dragSRPmoon.tikz'
            if not src.exists():
                print(f'  SKIP (missing): {src.name}')
                continue
            if PROCESS_TIKZ:
                fn = patch_anchor_tikz if 'anchors' in stem else preprocess_tikz
                title = fn(src, DST_RUNTIME / src.name)
            else:
                title = extract_title(src.read_text(encoding='utf-8'))
            caption = clean_caption(title)
            if not caption:
                metric, param = stem.split('_', 1)
                caption = (f"{METRIC_LABEL.get(metric, metric)} vs.\\ "
                           f"{PARAM_LABEL.get(param, param)}.")
            write_figure(f, 'Figures/runtime', src.name,
                         caption, f'fig:rt_{stem}')
            print(f'  runtime/{src.name}')

    # Alpha sweep tables
    f.write('\\subsubsection{Tables}\n\n')
    for fname, label in ALPHA_TABLES:
        src = SRC_TABLES / fname
        if not src.exists():
            print(f'  SKIP (missing): {fname}')
            continue
        if PROCESS_TIKZ:
            copy_table(src, DST_TABLES / fname)
        tbl_label = 'tbl:alpha_' + fname.replace('alpha_', '').replace('_J6dragSRPmoon.tex', '').replace('_', ':')
        write_table(f, fname, label, tbl_label)
        print(f'  table: {fname}')


# ---------------------------------------------------------------------------
# Section 3 — sensitivity sweep
# ---------------------------------------------------------------------------

SENS_ORBITS = ['LEO-1', 'MEO-1', 'HEO-E07', 'GEO-1']
SENS_METRICS = {
    'eff':    r'efficiency --- fevals vs.\ position error',
    'timing': r'timing --- wall time vs.\ position error',
}


def write_section_3(f):
    f.write('% ================================================================\n')
    f.write('\\subsection{Sensitivity Sweep}\n\n')

    # Sensitivity legend
    for orbit in SENS_ORBITS:
        src_leg = SRC_SENSITIVITY / f'{orbit}_J6dragSRPmoon_eff.tikz'
        if src_leg.exists():
            if PROCESS_TIKZ:
                build_legend(src_leg, DST_LEGEND_SENS, marker_scale=SENS_MARKER_SCALE)
            write_legend_figure(f, 'Figures', 'legend_sensitivity.tikz')
            break
    else:
        print('  WARNING: no eff tikz found for sensitivity legend.')

    sens_table_suffixes = [
        ('',        'Position error at $t_\\mathrm{end}$ (km)'),
        ('_wtime',  'Wall time (s)'),
        ('_fevals', 'Function evaluations'),
    ]

    for orbit in SENS_ORBITS:
        f.write(f'\\subsubsection{{{orbit}}}\n\n')
        for metric, metric_label in SENS_METRICS.items():
            src = SRC_SENSITIVITY / f'{orbit}_J6dragSRPmoon_{metric}.tikz'
            if not src.exists():
                print(f'  SKIP (missing): {src.name}')
                continue
            if PROCESS_TIKZ:
                preprocess_tikz(src, DST_SENSITIVITY / src.name,
                                marker_scale=SENS_MARKER_SCALE)
            caption = f'{orbit} --- {metric_label}.'
            write_figure(f, 'Figures/sensitivity', src.name,
                         caption, f'fig:sens_{slug(orbit)}_{metric}')
            print(f'  sensitivity/{src.name}')
        for suffix, tlabel in sens_table_suffixes:
            fname = f'sens_{orbit}_J6dragSRPmoon{suffix}.tex'
            src_t = SRC_TABLES / fname
            if not src_t.exists():
                print(f'  SKIP (missing): {fname}')
                continue
            if PROCESS_TIKZ:
                copy_table(src_t, DST_TABLES / fname)
            tbl_label = f'tbl:sens_{slug(orbit)}_{suffix.lstrip("_") or "r_err"}'
            write_table(f, fname, tlabel, tbl_label)
            print(f'  table: {fname}')


# ---------------------------------------------------------------------------
# Section 4 — full sweep (selected orbits in tex, all processed)
# ---------------------------------------------------------------------------

def write_section_4(f):
    f.write('% ================================================================\n')
    f.write('\\subsection{Full Sweep}\n\n')

    # Build legend from E07 (includes alpha, delta, N, and ode45/MPCM entries)
    src_leg = SRC_FULL / 'full_HEO-E07_J6dragSRPmoon.tikz'
    if src_leg.exists():
        if PROCESS_TIKZ:
            build_legend(src_leg, DST_LEGEND_FULL, label='full legend',
                         marker_scale=MARKER_SCALE)
        write_legend_figure(f, 'Figures', 'legend_full.tikz')
    else:
        print('  WARNING: full legend source not found.')

    # Process ALL full tikz files but only emit tex for FULL_TEX_ORBITS
    for src in sorted(SRC_FULL.glob('full_*.tikz')):
        orbit = src.stem.replace('full_', '').replace('_J6dragSRPmoon', '')
        if PROCESS_TIKZ:
            preprocess_tikz(src, DST_FULL / src.name, marker_scale=MARKER_SCALE)
            print(f'  full/{src.name}')
        if orbit in FULL_TEX_ORBITS:
            title = extract_title(src.read_text(encoding='utf-8'))
            caption = clean_caption(title) or f'Full sweep --- {orbit}.'
            write_figure(f, 'Figures/full', src.name,
                         caption, f'fig:full_{slug(orbit)}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if PROCESS_TIKZ:
        DST_ROOT.mkdir(parents=True, exist_ok=True)
        DST_TABLES.mkdir(parents=True, exist_ok=True)
        # Build alpha legend (MATLAB's legend.tikz has only \addlegendentry with no
        # \addplot, which breaks pgfplots at \end{axis})
        src_data = next(iter(SRC_ERROR_TIME.glob('*.tikz')), None)
        if src_data and SRC_LEGEND_ALPHA.exists():
            build_alpha_legend(src_data, SRC_LEGEND_ALPHA, DST_LEGEND_ALPHA)
        elif SRC_LEGEND_ALPHA.exists():
            shutil.copy2(SRC_LEGEND_ALPHA, DST_LEGEND_ALPHA)
            print(f'WARNING: no data tikz found; copied raw alpha legend -> {DST_LEGEND_ALPHA}')
        else:
            print(f'WARNING: alpha legend not found at {SRC_LEGEND_ALPHA}')

    with TEX_OUT.open('w', encoding='utf-8') as f:
        f.write('% Auto-generated by gen_figures_tex.py — do not edit by hand\n')
        f.write('% Run: python gen_figures_tex.py\n\n')
        write_section_1(f)
        write_section_2(f)
        write_section_3(f)
        write_section_4(f)

    print(f'\nDone.')
    print(f'  LaTeX fragment : {TEX_OUT}')
    if PROCESS_TIKZ:
        print(f'  Processed tikz : {DST_ROOT}')


if __name__ == '__main__':
    main()
