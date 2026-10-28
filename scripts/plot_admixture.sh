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
#      sample order: grouped by Super-population (then Population), labelled underneath.
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
for c in ("IID", "Population", "Super-population"):
    if c not in t.columns:
        sys.exit(f"ERROR: table lacks column {c!r}")
t = t.drop_duplicates("IID").set_index("IID")
meta = t.reindex(fam["IID"])[["Population", "Super-population"]].reset_index(drop=True)
missing = meta["Super-population"].isna().sum()
if missing:
    print(f"WARNING: {missing} fam samples are not in the table; grouped as 'unknown'")
meta = meta.fillna("unknown")

# Super-population order: roughly west to east / north to south; any others go at the end alphabetically.
ORDER = ["Iran", "Anatolia", "Caucuses", "Levant", "Arabia", "North-Africa", "Africa", "Southern-Europe",
         "Western-Europe", "Northern-Europe", "Eastern-Europe", "Central-Asia", "South-Asia",
         "East-Asia", "South-East-Asia", "Oceania", "Americas", "unknown"]
rank = {s: i for i, s in enumerate(ORDER)}
meta["rank"] = meta["Super-population"].map(lambda s: rank.get(s, len(ORDER) - 1))
order = meta.assign(i=np.arange(n)).sort_values(["rank", "Super-population", "Population", "i"])["i"].to_numpy()
groups = meta["Super-population"].to_numpy()[order]
bounds = [0] + [i for i in range(1, n) if groups[i] != groups[i - 1]] + [n]

# Component colours are by index only (ADMIXTURE components are not matched between Ks).
PALETTE = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffd92f", "#a65628", "#f781bf",
           "#17becf", "#999999", "#1b9e77", "#7570b3", "#66a61e", "#e6ab02", "#666666", "#8dd3c7"]

fig, axes = plt.subplots(len(ks), 1, figsize=(18, 1.1 * len(ks) + 1.6), sharex=True)
axes = np.atleast_1d(axes)
edges = np.arange(n + 1) - 0.5
for ax, k in zip(axes, ks):
    q = Q[k][order]
    # step-stacked areas (one pixel column per sample); much faster than thousands of bar patches
    ax.stackplot(edges, np.vstack([q, q[-1:]]).T, step="post",
                 colors=[PALETTE[j % len(PALETTE)] for j in range(k)], linewidth=0)
    for b in bounds[1:-1]:
        ax.axvline(b - 0.5, color="black", lw=0.8)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel(f"K={k}", rotation=0, ha="right", va="center", fontsize=10)
    for s in ax.spines.values():
        s.set_linewidth(0.6)
mids = [(a + b) / 2 - 0.5 for a, b in zip(bounds[:-1], bounds[1:])]
axes[-1].set_xticks(mids)
axes[-1].set_xticklabels([f"{groups[a]} (n={b - a})" for a, b in zip(bounds[:-1], bounds[1:])],
                         rotation=60, ha="right", fontsize=9)
axes[-1].tick_params(axis="x", length=0)
fig.subplots_adjust(left=0.05, right=0.995, top=0.99, bottom=0.13, hspace=0.12)
for ext in ("png", "pdf"):
    fig.savefig(f"{prefix}.admixture_stacked.{ext}", dpi=200)
plt.close(fig)
print(f"stacked plot: K={ks[0]}..{ks[-1]}, {n} samples, {len(bounds) - 1} super-population blocks")
PY

echo "done"
