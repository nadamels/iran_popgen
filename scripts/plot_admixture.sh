#!/bin/bash

#SBATCH -c 2

# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/glimpse2

# Plot ADMIXTURE results:
#   1. CV error against K (from the "CV error (K=..): ..." lines of the ADMIXTURE logs;
#      the K with the lowest CV error is marked).
#   2. One stacked-bar ("structure") plot with every K, one panel per K, all panels using the same
#      sample order. Left to right: ancient samples grouped by archaeological period (oldest first;
#      ancient Iranians last), a THICK black line, then modern samples (modern Iranians first, then the
#      rest by Super-population). Ancient = has an Archeological Date or Period in the table.
#
# Usage:   sbatch plot_admixture.sh <prefix> <cv_logs> [mega_sample_table.xlsx|.csv]
#          bash   plot_admixture.sh <prefix> <cv_logs> [mega_sample_table.xlsx|.csv]
#
#   <prefix>   expects  <prefix>.fam  and  <prefix>.<K>.Q  for each K found
#              (e.g. q_files_and_fam/all_anc_mod.run3.LDpruned)
#   <cv_logs>  ADMIXTURE stdout log file(s): a file, a directory (all *.log / *.out inside) or a glob
#              in quotes. Use "none" to skip the CV plot. If a K appears in several logs (e.g. several
#              seeds) the CV errors are shown as points with the mean joined by a line.
#   table      .xlsx (sheet "mega_table") or .csv. Default: $MASTER_TABLE, else
#              metadata/mega_sample_table.xlsx under $SLURM_SUBMIT_DIR (or the current directory),
#              falling back to the .csv next to it if openpyxl is missing
#
# Outputs (next to the prefix): <prefix>.cv_error.{png,pdf}, <prefix>.cv_error.tsv,
#                               <prefix>.admixture_stacked.{png,pdf}

set -euo pipefail

PREFIX="${ARGS[0]:-}"
CV_LOGS="${ARGS[1]:-}"
if [[ -z "$PREFIX" || -z "$CV_LOGS" ]]; then
    echo "Usage: $0 <prefix> <cv_logs|none> [mega_sample_table.xlsx|.csv]" >&2
    exit 1
fi

BASE="${SLURM_SUBMIT_DIR:-$PWD}"
TABLE="${ARGS[2]:-${MASTER_TABLE:-$BASE/metadata/mega_sample_table.xlsx}}"

if [[ "$TABLE" == *.xlsx || "$TABLE" == *.xls ]] && ! python3 -c "import openpyxl" 2>/dev/null; then
    alt="${TABLE%.*}.csv"
    if [[ -s "$alt" ]]; then
        echo "NOTE: openpyxl not available; using $alt instead of $TABLE" >&2
        TABLE="$alt"
    fi
fi

for f in "${PREFIX}.fam" "$TABLE"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done

python3 -c "import numpy, pandas, matplotlib" 2>/dev/null || {
    echo "ERROR: python3 in the active environment needs numpy, pandas and matplotlib" >&2; exit 1; }

echo "prefix : $PREFIX"
echo "cv logs: $CV_LOGS"
echo "table  : $TABLE"

python3 - "$PREFIX" "$CV_LOGS" "$TABLE" <<'PY'
import glob
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

prefix, cv_logs, table_path = sys.argv[1], sys.argv[2], sys.argv[3]

# ---- 1. CV error vs K ---------------------------------------------------------------------------
CV = re.compile(r"CV error \(K=(\d+)\):\s*([0-9.eE+-]+)")
if cv_logs.lower() != "none":
    if os.path.isdir(cv_logs):
        files = sorted(glob.glob(os.path.join(cv_logs, "*.log")) + glob.glob(os.path.join(cv_logs, "*.out")))
    else:
        files = sorted(glob.glob(cv_logs))
    rows = []
    for fn in files:
        with open(fn, errors="replace") as fh:
            rows += [(int(k), float(v), os.path.basename(fn)) for k, v in CV.findall(fh.read())]
    if not rows:
        sys.exit(f"ERROR: no 'CV error (K=..)' lines found in {len(files)} file(s) matching {cv_logs!r}")
    cv = pd.DataFrame(rows, columns=["K", "CV_error", "log"]).sort_values(["K", "log"])
    cv.to_csv(f"{prefix}.cv_error.tsv", sep="\t", index=False)
    mean = cv.groupby("K")["CV_error"].mean()
    best = int(mean.idxmin())
    print(f"CV error: {len(cv)} values from {len(files)} log(s), K={mean.index.min()}..{mean.index.max()}; "
          f"lowest mean CV error at K={best} ({mean[best]:.5f})")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(mean.index, mean.values, "-", color="0.4", lw=1.5, zorder=1)
    ax.scatter(cv["K"], cv["CV_error"], s=28, color="#1f4e79", zorder=2)
    ax.scatter([best], [mean[best]], s=140, facecolor="none", edgecolor="#c00000", lw=2, zorder=3,
               label=f"lowest: K={best}")
    ax.set_xticks(mean.index)
    ax.set_xlabel("K")
    ax.set_ylabel("Cross-validation error")
    ax.legend(frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{prefix}.cv_error.{ext}", dpi=300)
    plt.close(fig)

# ---- 2. stacked plot of all K -------------------------------------------------------------------
qfiles = {}
for fn in glob.glob(f"{glob.escape(prefix)}.*.Q"):
    m = re.fullmatch(r"(\d+)", os.path.basename(fn)[len(os.path.basename(prefix)) + 1:-2])
    if m:
        qfiles[int(m.group(1))] = fn
if not qfiles:
    sys.exit(f"ERROR: no {prefix}.<K>.Q files found")
ks = sorted(qfiles)

fam = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None, dtype=str,
                  names=["FID", "IID", "PID", "MID", "SEX", "PHENO"])
n = len(fam)
Q = {}
for k in ks:
    q = np.loadtxt(qfiles[k], ndmin=2)
    if q.shape != (n, k):
        sys.exit(f"ERROR: {qfiles[k]} has shape {q.shape}, expected ({n}, {k}) from the fam")
    Q[k] = q

if table_path.lower().endswith((".xlsx", ".xls")):
    try:
        t = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
    except ImportError:
        sys.exit("ERROR: reading .xlsx needs openpyxl; pass the .csv instead")
else:
    t = pd.read_csv(table_path, dtype=str)
for c in ("IID", "Population", "Super-population", "Archeological Period", "Archeological Date"):
    if c not in t.columns:
        sys.exit(f"ERROR: table lacks column {c!r}")
t = t.drop_duplicates("IID").set_index("IID")
cols = ["Population", "Super-population", "Archeological Period", "Archeological Date"]
meta = t.reindex(fam["IID"])[cols].reset_index(drop=True)
missing = meta["Super-population"].isna().sum()
if missing:
    print(f"WARNING: {missing} fam samples are not in the table; grouped as 'unknown'")
meta["Super-population"] = meta["Super-population"].fillna("unknown")

# Ancient = has an archaeological date or period. Layout, left to right:
#   ancient non-Iranians by period (oldest first), ancient Iranians by period (at the end of the ancient
#   block) | THICK LINE | modern Iranians first, then the other modern samples by super-population.
meta["ancient"] = meta["Archeological Date"].notna() | meta["Archeological Period"].notna()
meta["period"] = meta["Archeological Period"].fillna("unknown")
PERIODS = ["PPN", "Natufian", "M / N", "N", "EN", "C / EN", "ChL", "EMBA", "BA", "MLBA", "LBA", "LBA / IA",
           "IA", "Urartian", "Rom", "RomByz", "Byz", "Mdv", "Anc", "unknown"]
prank = {s: i for i, s in enumerate(PERIODS)}
unk = [p for p in meta.loc[meta["ancient"], "period"].unique() if p not in prank]
if unk:
    print(f"WARNING: periods not in the ordering list, placed last: {sorted(unk)}")
meta["prank"] = meta["period"].map(lambda s: prank.get(s, len(PERIODS)))
meta["iran"] = meta["Super-population"].eq("Iran")

# Super-population order for moderns: roughly west to east / north to south; others at the end.
ORDER = ["Iran", "Anatolia", "Caucuses", "Levant", "Arabia", "North-Africa", "Africa", "Southern-Europe",
         "Western-Europe", "Northern-Europe", "Eastern-Europe", "Central-Asia", "South-Asia",
         "East-Asia", "South-East-Asia", "Oceania", "Americas", "unknown"]
rank = {s: i for i, s in enumerate(ORDER)}
meta["rank"] = meta["Super-population"].map(lambda s: rank.get(s, len(ORDER) - 1))

anc = meta["ancient"]
# block label and sort keys: (ancient first, ancient-Iran last within ancients, then period / super-pop)
meta["label"] = np.where(anc, meta["Super-population"] + " · " + meta["period"],
                         meta["Super-population"].where(~meta["iran"], "Iran (modern)"))
meta["k1"] = np.where(anc, 0, 1)
meta["k2"] = np.where(anc, meta["iran"].astype(int), 0)
meta["k3"] = np.where(anc, meta["prank"], meta["rank"])
meta["i"] = np.arange(n)
order = meta.sort_values(["k1", "k2", "k3", "Super-population", "Population", "i"])["i"].to_numpy()
groups = meta["label"].to_numpy()[order]
bounds = [0] + [i for i in range(1, n) if groups[i] != groups[i - 1]] + [n]
n_anc = int(anc.sum())
print(f"ancient: {n_anc} samples; modern: {n - n_anc} samples "
      f"({int((~anc & meta['iran']).sum())} modern Iran-labelled placed first among moderns)")

# Component colours are by index only (ADMIXTURE components are not matched between Ks).
PALETTE = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffd92f", "#a65628", "#f781bf",
           "#17becf", "#999999", "#1b9e77", "#7570b3", "#66a61e", "#e6ab02", "#666666", "#8dd3c7"]

# Two panels per K: ancient (left) and modern (right), joined by a thick black line. The ancient panel is
# drawn wider than its sample share (ANC_SCALE x) so its many small period blocks stay readable.
n_mod = n - n_anc
has_anc = 0 < n_anc < n
FIG_W = 18
LABEL_PT = 7
if has_anc:
    n_anc_blocks = sum(1 for a in bounds[:-1] if a < n_anc)
    need_in = n_anc_blocks * LABEL_PT * 1.6 / 72                  # room for one rotated label per block
    frac = min(max(need_in / (FIG_W * 0.945), 0.2), 0.75)          # share of the width for ancients
    ANC_SCALE = frac / (1 - frac) * n_mod / n_anc
    print(f"ancient panel drawn {ANC_SCALE:.1f}x wider than its sample share ({n_anc_blocks} blocks)")
fig, axes = plt.subplots(len(ks), 2, figsize=(FIG_W, 1.1 * len(ks) + 2.4),
                         gridspec_kw={"width_ratios": [n_anc * ANC_SCALE, n_mod] if has_anc else [1, 1e-9],
                                      "wspace": 0})
edges = np.arange(n + 1) - 0.5
spans = [(0, n_anc), (n_anc, n)]
for r, k in enumerate(ks):
    q = Q[k][order]
    for ax, (lo, hi) in zip(axes[r], spans):
        if hi <= lo:
            ax.axis("off")
            continue
        # step-stacked areas (one pixel column per sample); much faster than thousands of bar patches
        ax.stackplot(edges[lo:hi + 1], np.vstack([q[lo:hi], q[hi - 1:hi]]).T, step="post",
                     colors=[PALETTE[j % len(PALETTE)] for j in range(k)], linewidth=0)
        for b in bounds[1:-1]:
            if lo < b < hi:
                ax.axvline(b - 0.5, color="black", lw=0.8)
        ax.set_xlim(lo - 0.5, hi - 0.5)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.6)
    axes[r][0].set_ylabel(f"K={k}", rotation=0, ha="right", va="center", fontsize=10)
    if has_anc:                                                  # thick ancient | modern line
        axes[r][0].spines["right"].set_linewidth(3.5)
        axes[r][1].spines["left"].set_linewidth(3.5)

def spread(pos, gap, lo, hi):
    """Nudge sorted label positions apart to at least `gap`, keeping them inside [lo, hi]."""
    pos = list(pos)
    for i in range(1, len(pos)):
        pos[i] = max(pos[i], pos[i - 1] + gap)
    pos[-1] = min(pos[-1], hi)
    for i in range(len(pos) - 2, -1, -1):
        pos[i] = min(pos[i], pos[i + 1] - gap)
    pos[0] = max(pos[0], lo)
    return pos

for ax, (lo, hi) in zip(axes[-1], spans):
    if hi <= lo:
        continue
    blocks = [(a, b) for a, b in zip(bounds[:-1], bounds[1:]) if lo <= a and b <= hi]
    mids = [(a + b) / 2 - 0.5 for a, b in blocks]
    width_in = ax.get_position().width * FIG_W
    gap = (hi - lo) * (LABEL_PT * 1.3 / 72) / width_in           # one text line height, in data units
    lpos = spread(mids, gap, lo, hi - 1)
    for (a, b), m, lx in zip(blocks, mids, lpos):
        ax.annotate(f"{groups[a]} (n={b - a})", xy=(m, 0), xycoords=("data", "axes fraction"),
                    xytext=(lx, -0.06), textcoords=("data", "axes fraction"), rotation=90,
                    ha="center", va="top", fontsize=LABEL_PT,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="0.4", shrinkA=0, shrinkB=0))
if has_anc:
    axes[0][0].text(0.5, 1.04, "Ancient", ha="center", va="bottom", fontsize=11, fontweight="bold",
                    transform=axes[0][0].transAxes)
    axes[0][1].text(0.5, 1.04, "Modern", ha="center", va="bottom", fontsize=11, fontweight="bold",
                    transform=axes[0][1].transAxes)
fig.subplots_adjust(left=0.05, right=0.995, top=0.96, bottom=0.17, hspace=0.12)
for ext in ("png", "pdf"):
    fig.savefig(f"{prefix}.admixture_stacked.{ext}", dpi=200)
plt.close(fig)
print(f"stacked plot: K={ks[0]}..{ks[-1]}, {n} samples, {len(bounds) - 1} blocks")
PY

echo "done"
