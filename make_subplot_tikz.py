#!/usr/bin/env python3
"""
make_subplot_tikz.py

Further-process already-processed full-sweep tikz files for a 2x2 subplot panel.

Source: results/figures_processed/full/full_{ORBIT}_J6dragSRPmoon.tikz
Output: results/figures_processed/full/subplot_{ORBIT}_J6dragSRPmoon.tikz

Transformations applied:
  - width  -> 0.45\\linewidth
  - height -> 0.3\\linewidth
  - mark size values halved (x0.5)
  - ode45 \\addplot blocks commented out
  - "MPCM failed" -> "MGPI failed" in legend entries

Also generates legend_full_sub.tikz from legend_full.tikz:
  - mark sizes halved
  - legend columns=4
  - wrapped in \\fbox

Usage:
    python make_subplot_tikz.py E03 E05 E07 E09
"""

import re
import sys
from pathlib import Path

from gen_figures_tex import scale_marker_sizes

# ---------------------------------------------------------------------------
SUBPLOT_WIDTH  = r'0.45\linewidth'
SUBPLOT_HEIGHT = r'0.3\linewidth'
MARKER_SCALE   = 0.5
FORCE_MODEL    = 'J6dragSRPmoon'

SRC_FULL   = Path('results/figures_processed/full')
DST_FULL   = Path('results/figures_processed/full')
LEGEND_SRC = Path('results/figures_processed/legend_full.tikz')
LEGEND_DST = Path('results/figures_processed/legend_full_sub.tikz')
# ---------------------------------------------------------------------------


def _expand_orbit(name: str) -> str:
    """Accept E03 shorthand or full HEO-E03 form."""
    for prefix in ('HEO-', 'LEO-', 'MEO-', 'GEO-'):
        if name.startswith(prefix):
            return name
    candidate = 'HEO-' + name
    if (SRC_FULL / f'full_{candidate}_{FORCE_MODEL}.tikz').exists():
        return candidate
    return name


def comment_ode45_blocks(content: str) -> str:
    """Comment \\addplot blocks whose associated \\addlegendentry (active or
    already commented) contains 'ode45'.
    Rename 'MPCM failed' -> 'MGPI failed' in both active and commented legend entries."""
    lines = content.splitlines(keepends=True)
    result = []
    i = 0
    while i < len(lines):
        line = lines[i]

        # Rename MPCM -> MGPI (active or already-commented legend entry)
        if re.match(r'\s*%?\\addlegendentry', line) and 'MPCM' in line:
            result.append(line.replace('MPCM', 'MGPI'))
            i += 1
            continue

        # Detect start of a non-forget-plot \addplot block
        if re.match(r'\s*\\addplot', line) and 'forget plot' not in line:
            block = [line]
            # Single-line block ends with }; on the same line
            if line.rstrip().endswith('};'):
                j = i + 1
            else:
                j = i + 1
                while j < len(lines):
                    block.append(lines[j])
                    if re.match(r'\s*\};', lines[j].strip()):
                        j += 1
                        break
                    j += 1
            # Peek past blank lines for legend entry (active OR commented)
            k = j
            while k < len(lines) and lines[k].strip() == '':
                k += 1
            if k < len(lines) and re.match(r'\s*%?\\addlegendentry', lines[k]) and 'ode45' in lines[k]:
                result.extend('%' + bl for bl in block)
                result.extend(lines[j:k])
                result.append(lines[k] if lines[k].lstrip().startswith('%') else '%' + lines[k])
                i = k + 1
            else:
                result.extend(block)
                i = j
            continue

        result.append(line)
        i += 1
    return ''.join(result)


def process_subplot(src: Path, dst: Path):
    content = src.read_text(encoding='utf-8')

    # 1. Comment ode45 blocks, rename MPCM->MGPI
    content = comment_ode45_blocks(content)

    # 2. Replace width and height
    content = re.sub(r'(\s*)width\s*=\s*[^,\n]+,',
                     lambda m: f'{m.group(1)}width={SUBPLOT_WIDTH},', content)
    content = re.sub(r'(\s*)height\s*=\s*[^,\n]+,',
                     lambda m: f'{m.group(1)}height={SUBPLOT_HEIGHT},', content)

    # 3. Halve marker sizes
    content = scale_marker_sizes(content, MARKER_SCALE)

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(content, encoding='utf-8')


def make_subplot_legend():
    """Generate legend_full_sub.tikz: halved markers, 4 columns, fbox, no ode45 entries."""
    if not LEGEND_SRC.exists():
        print(f'  SKIP legend (not found): {LEGEND_SRC}')
        return
    content = LEGEND_SRC.read_text(encoding='utf-8')
    # Comment ode45 addplot+addlegendentry pairs (active in legend file)
    content = comment_ode45_blocks(content)
    content = scale_marker_sizes(content, MARKER_SCALE)
    content = content.replace('legend columns=-1', 'legend columns=4')
    content = content.replace(r'\begin{tikzpicture}',
                              r'\fbox{' + '\n' + r'\begin{tikzpicture}')
    content = content.replace(r'\end{tikzpicture}',
                              r'\end{tikzpicture}' + '\n' + '}')
    LEGEND_DST.parent.mkdir(parents=True, exist_ok=True)
    LEGEND_DST.write_text(content, encoding='utf-8')
    print(f'  legend -> {LEGEND_DST.name}')


def main():
    if len(sys.argv) < 2:
        print('Usage: python make_subplot_tikz.py ORBIT [ORBIT ...]')
        print('  e.g. python make_subplot_tikz.py E03 E05 E07 E09')
        sys.exit(1)

    for arg in sys.argv[1:]:
        orbit = _expand_orbit(arg)
        src = SRC_FULL / f'full_{orbit}_{FORCE_MODEL}.tikz'
        dst = DST_FULL / f'subplot_{orbit}_{FORCE_MODEL}.tikz'
        if not src.exists():
            print(f'  SKIP (not found): {src}')
            continue
        process_subplot(src, dst)
        print(f'  {orbit} -> {dst.name}')

    make_subplot_legend()


if __name__ == '__main__':
    main()
