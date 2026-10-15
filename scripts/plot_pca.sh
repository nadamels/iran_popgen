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

# Plot a PCA from the mega sample table.
#   COLOUR = Super-population (colours from Plotting/pca_colors.txt)
#   SIZE   = age of the sample (older = bigger), ancient samples only
#   Modern samples = tiny dots (coloured by super-population); no marker shapes for now.
#   Modern IRANIANS (Super-population Iran, the target population) = large outlined SQUARES drawn
#   on top of everything else.
#
# Usage:   sbatch plot_pca.sh <prefix> [mega_sample_table.xlsx|.csv] [pca_colors.txt] [PC1_range] [PC2_range]
#          bash   plot_pca.sh <prefix> [mega_sample_table.xlsx|.csv] [pca_colors.txt] [PC1_range] [PC2_range]
#
#   Zoom example (PC1 from -0.01 to 0.01, PC2 from 0.0 to 0.02):
#          sbatch plot_pca.sh <prefix> <table> <colors> -0.01,0.01 0.0,0.02
#   A range is "min,max" (or "min:max"); use "auto" (or "-") to leave that axis unzoomed.
#   Defaults for the first three arguments can be skipped with "" (empty string).
#
#   <prefix>  expects  <prefix>.fam  <prefix>.eigenvec  <prefix>.eigenval
#             (PLINK 1.9 or PLINK 2 eigenvec; with or without a header line)
#   table     .xlsx (sheet "mega_table"; needs openpyxl) or .csv. Default: $MASTER_TABLE, else
#             metadata/mega_sample_table.xlsx under $SLURM_SUBMIT_DIR (or the current directory),
#             falling back to the .csv next to it if openpyxl is missing
#   colours   "name: #hex" lines. Default: $COLORS_FILE, else Plotting/pca_colors.txt (same base dir)
#
# Age: parsed from 'Archeological Date' (e.g. "5898-5531 calBCE", "600-800 CE", "1500-800 BCE")
#      as the range midpoint in years before present (1950). Unparseable dates are drawn at the
#      smallest ancient size and counted in the log.
#
#   When a PC1/PC2 range is given, ONLY the zoomed PC1 vs PC2 plot is made, written to
#   <prefix>.pca_PC1_PC2.zoom_x<min>_<max>_y<min>_<max>.{png,pdf} (the full plots are not overwritten).
#
# Outputs (next to the prefix):  <prefix>.pca_PC1_PC2.{png,pdf}, <prefix>.pca_PC3_PC4.{png,pdf},
#                                <prefix>.pca_samples.tsv (IID, super-population, age BP, colour)

set -euo pipefail

PREFIX="${ARGS[0]:-}"
if [[ -z "$PREFIX" ]]; then
    echo "Usage: $0 <prefix> [mega_sample_table.xlsx|.csv] [pca_colors.txt] [PC1_min,max] [PC2_min,max]" >&2
    exit 1
fi
BASE="${SLURM_SUBMIT_DIR:-$PWD}"
TABLE="${ARGS[1]:-${MASTER_TABLE:-$BASE/metadata/mega_sample_table.xlsx}}"
COLORS="${ARGS[2]:-${COLORS_FILE:-$BASE/Plotting/pca_colors.txt}}"
PC1_RANGE="${ARGS[3]:-auto}"
PC2_RANGE="${ARGS[4]:-auto}"

if [[ "$TABLE" == *.xlsx || "$TABLE" == *.xls ]] && ! python3 -c "import openpyxl" 2>/dev/null; then
    alt="${TABLE%.*}.csv"
    if [[ -s "$alt" ]]; then
        echo "NOTE: openpyxl not available; using $alt instead of $TABLE" >&2
        TABLE="$alt"
    fi
fi

for f in "${PREFIX}.fam" "${PREFIX}.eigenvec" "${PREFIX}.eigenval" "$TABLE" "$COLORS"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done

python3 -c "import numpy, pandas, matplotlib" 2>/dev/null || {
    echo "ERROR: python3 in the active environment needs numpy, pandas and matplotlib" >&2; exit 1; }

echo "prefix : $PREFIX"
echo "table  : $TABLE"
echo "colours: $COLORS"
echo "PC1 range: $PC1_RANGE   PC2 range: $PC2_RANGE"

python3 - "$PREFIX" "$TABLE" "$COLORS" "$PC1_RANGE" "$PC2_RANGE" <<'PY'
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

prefix, table_path, colors_path = sys.argv[1], sys.argv[2], sys.argv[3]

def parse_range(txt, name):
    """'min,max' or 'min:max' -> (min, max); 'auto' / '-' / '' -> None"""
    if txt.strip().lower() in ("", "auto", "-", "none"):
        return None
    parts = [x for x in re.split(r"[,:]", txt.strip()) if x != ""]
    try:
        lo, hi = (float(x) for x in parts)
    except ValueError:
        sys.exit(f"ERROR: {name} range {txt!r} must look like min,max (e.g. -0.01,0.01) or 'auto'")
    if not lo < hi:
        sys.exit(f"ERROR: {name} range {txt!r}: min must be smaller than max")
    return lo, hi
xr, yr = parse_range(sys.argv[4], "PC1"), parse_range(sys.argv[5], "PC2")
zoom = xr is not None or yr is not None

# ---- read fam / eigenvec / eigenval ------------------------------------------------------------
fam = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None, dtype=str,
                  names=["FID", "IID", "PID", "MID", "SEX", "PHENO"])
# Un-repaired fam where every IID is one constant string (ID stored in FID): fall back to FID.
use_fid = bool(fam["IID"].duplicated().any() and not fam["FID"].duplicated().any())
if use_fid:
    print("WARNING: IIDs in the fam are not unique but FIDs are; using FID as the sample ID")
fam["sid"] = fam["FID"] if use_fid else fam["IID"]

with open(f"{prefix}.eigenvec") as fh:
    first = fh.readline().split()
has_header = first and first[0].startswith("#")
if has_header:                                   # PLINK 2: "#FID IID PC1..." or "#IID PC1..."
    ev = pd.read_csv(f"{prefix}.eigenvec", sep=r"\s+", dtype=str)
    ev.columns = [c.lstrip("#") for c in ev.columns]
    pc_cols = [c for c in ev.columns if c.upper().startswith("PC")]
    id_cols = [c for c in ev.columns if c not in pc_cols]
    ev["FID"] = ev["FID"] if "FID" in ev else "0"
else:                                            # PLINK 1.9: FID IID PC1 PC2 ...
    ev = pd.read_csv(f"{prefix}.eigenvec", sep=r"\s+", header=None, dtype=str)
    pc_cols = [f"PC{i}" for i in range(1, ev.shape[1] - 1)]
    ev.columns = ["FID", "IID"] + pc_cols
ev[pc_cols] = ev[pc_cols].astype(float)

key = "FID" if use_fid else "IID"
ev = ev.rename(columns={key: "sid"})
if len(pc_cols) < 2:
    sys.exit("ERROR: eigenvec has fewer than 2 PCs")
eigval = np.loadtxt(f"{prefix}.eigenval", ndmin=1)
pct = 100 * eigval / eigval.sum()

# ---- fam <-> eigenvec consistency --------------------------------------------------------------
fam_ids, ev_ids = set(fam["sid"]), set(ev["sid"])
if fam_ids != ev_ids:
    print(f"WARNING: fam has {len(fam_ids)} samples, eigenvec {len(ev_ids)}; "
          f"{len(fam_ids - ev_ids)} only in fam, {len(ev_ids - fam_ids)} only in eigenvec")
if ev["sid"].duplicated().any():
    sys.exit("ERROR: duplicate sample IDs in eigenvec; cannot join to the sample table unambiguously")
ev = ev[ev["sid"].isin(fam_ids)]               # fam defines the sample set

# ---- mega table: super-population, age ---------------------------------------------------------
if table_path.lower().endswith((".xlsx", ".xls")):
    try:
        t = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
    except ImportError:
        sys.exit("ERROR: reading .xlsx needs openpyxl; pass the .csv instead")
else:
    t = pd.read_csv(table_path, dtype=str)
need = ["IID", "Super-population", "Archeological Date", "Undergone Imputation"]
miss_cols = [c for c in need if c not in t.columns]
if miss_cols:
    sys.exit(f"ERROR: table lacks column(s) {miss_cols}")
t = t.drop_duplicates("IID")

DATE = re.compile(r"^\s*(?:(?P<a>\d+)\s*-\s*(?P<b>\d+)|(?P<c>\d+))\s*(?P<era>cal\s*BCE|BCE|cal\s*CE|CE|AD|cal\s*AD|calBP|BP)", re.I)
def age_bp(txt):
    """Midpoint of the date range in years before present (1950)."""
    m = DATE.match(txt) if isinstance(txt, str) else None
    if not m:
        return np.nan
    lo, hi = (m["a"], m["b"]) if m["a"] else (m["c"], m["c"])
    mid = (float(lo) + float(hi)) / 2
    era = re.sub(r"\s+", "", m["era"]).lower()
    if era.endswith("bce"):
        return mid + 1950
    if era.endswith("bp"):
        return mid
    return 1950 - mid                                   # CE / AD

t["is_ancient"] = t["Undergone Imputation"].str.strip().str.lower() == "yes"
t["age"] = [age_bp(d) if a else np.nan for d, a in zip(t["Archeological Date"], t["is_ancient"])]

df = ev.merge(t[["IID", "Super-population", "is_ancient", "age"]], left_on="sid", right_on="IID", how="left")
unk = df["IID"].isna()
if unk.any():
    print(f"WARNING: {unk.sum()} samples not in the mega table (plotted grey 'x'): "
          f"{', '.join(df.loc[unk, 'sid'].head(10))}{' ...' if unk.sum() > 10 else ''}")
df["is_ancient"] = df["is_ancient"].fillna(False).astype(bool)
df["sp"] = df["Super-population"].fillna("Unassigned")
n_bad = int((df["is_ancient"] & df["age"].isna()).sum())
if n_bad:
    print(f"NOTE: {n_bad} ancient samples have no parseable date; drawn at the smallest ancient size")

# ---- colours -----------------------------------------------------------------------------------
def norm(s):
    return re.sub(r"[^a-z]", "", s.lower())
file_cols = {}
with open(colors_path) as fh:
    for line in fh:
        if ":" in line:
            k, v = line.split(":", 1)
            if v.strip().startswith("#"):
                file_cols[norm(k)] = v.strip()
# table super-population name -> key in the colour file (names differ slightly)
ALIAS = {"iran": "iraniansamplesofthisstudy", "arabia": "arabian", "levant": "levantine",
         "northafrica": "africa"}
# not listed in pca_colors.txt: fallback colours (edit the file to override; add a line "Name: #hex")
FALLBACK = {"anatolia": "#8c564b", "americas": "#2ca02c", "oceania": "#17becf", "unassigned": "#7f7f7f"}
def colour_of(sp):
    k = norm(sp)
    for key in (k, ALIAS.get(k)):
        if key and key in file_cols:
            return file_cols[key], True
    return FALLBACK.get(k, "#7f7f7f"), False
sps = sorted(df["sp"].unique())
col_of, fallback_used = {}, []
for sp in sps:
    col_of[sp], ok = colour_of(sp)
    if not ok:
        fallback_used.append(sp)
if fallback_used:
    print("NOTE: no colour in pca_colors.txt for: " + ", ".join(fallback_used) +
          " (fallback colours used; add 'Name: #hex' lines to set them)")

# ---- size from age -----------------------------------------------------------------------------
S_MIN, S_MAX = 22.0, 300.0                      # marker area (pt^2) for the youngest / oldest sample
AGE_MAX = float(np.nanmax(df.loc[df.is_ancient, "age"])) if df.is_ancient.any() else 1.0
AGE_MAX = max(AGE_MAX, 1.0)
def size_of(age):
    a = np.where(np.isfinite(age), np.clip(age, 0, AGE_MAX), 0.0)
    return S_MIN + (S_MAX - S_MIN) * a / AGE_MAX
df["size"] = size_of(df["age"].to_numpy(dtype=float))

# ---- plot --------------------------------------------------------------------------------------
IRAN_SIZE = 110.0                               # marker area (pt^2) for modern Iranians (target population)
anc = df[df.is_ancient]
mod_all = df[~df.is_ancient & df["IID"].notna()]
is_iran = mod_all["sp"].str.lower() == "iran"
mod, mod_iran = mod_all[~is_iran], mod_all[is_iran]
def plot_pair(a, b, out, xlim=None, ylim=None):
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_axes([0.06, 0.07, 0.66, 0.88])
    # modern: tiny dots, bottom layer
    for sp, s in mod.groupby("sp"):
        ax.scatter(s[a], s[b], s=3, c=[col_of[sp]], marker="o", edgecolors="none", alpha=0.7, zorder=1)
    u = df[df["IID"].isna()]
    if len(u):
        ax.scatter(u[a], u[b], s=14, c="grey", marker="x", zorder=2)
    # ancient: colour = super-population, size = age; younger on top so big old points stay visible
    s = anc.sort_values("age", ascending=False, na_position="last")
    ax.scatter(s[a], s[b], s=s["size"], c=[col_of[x] for x in s["sp"]], marker="o",
               edgecolors="black", linewidths=0.5, alpha=0.9, zorder=3)
    # modern Iranians: large, on top of ancient points too
    if len(mod_iran):
        ax.scatter(mod_iran[a], mod_iran[b], s=IRAN_SIZE, c=[col_of[x] for x in mod_iran["sp"]], marker="s",
                   edgecolors="black", linewidths=0.8, alpha=0.95, zorder=5)
    ia, ib = pc_cols.index(a), pc_cols.index(b)
    ax.set_xlabel(f"{a} ({pct[ia]:.2f}% of listed eigenvalues)" if ia < len(pct) else a)
    ax.set_ylabel(f"{b} ({pct[ib]:.2f}% of listed eigenvalues)" if ib < len(pct) else b)
    title = f"{prefix.split('/')[-1]}  |  n={len(df)}  (ancient {len(anc)}, modern {len(mod) + len(mod_iran)}, of which Iranian {len(mod_iran)})"
    if xlim is not None or ylim is not None:
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        sel = np.ones(len(df), bool)
        if xlim is not None:
            sel &= ((df[a] >= xlim[0]) & (df[a] <= xlim[1])).to_numpy()
        if ylim is not None:
            sel &= ((df[b] >= ylim[0]) & (df[b] <= ylim[1])).to_numpy()
        title += f"\nZOOM: {a} {xlim if xlim is not None else 'auto'}, {b} {ylim if ylim is not None else 'auto'}  -  {int(sel.sum())} samples in view"
        print(f"samples in the zoom window: {int(sel.sum())} of {len(df)}")
    ax.set_title(title)
    ax.grid(alpha=0.2)

    # legend 1: colour = super-population (ancient / modern counts)
    h1 = []
    for sp in sorted(sps, key=lambda x: (x == "Unassigned", x)):
        na, nm = int(((anc.sp == sp)).sum()), int(((mod_all.sp == sp)).sum())
        h1.append(Line2D([], [], marker="o", ls="", color=col_of[sp], markersize=8,
                         markeredgecolor="black", markeredgewidth=0.4,
                         label=f"{sp} (ancient {na}, modern {nm})"))
    # legend 2: size = age
    ticks = [t_ for t_ in (500, 2000, 4000, 6000, 8000, 10000, 12000) if t_ <= AGE_MAX * 1.05]
    h2 = [Line2D([], [], marker="o", ls="", color="0.75", markeredgecolor="black", markeredgewidth=0.5,
                 markersize=float(np.sqrt(size_of(np.array([t_]))[0])), label=f"{t_:,} years BP")
          for t_ in ticks]
    h2.append(Line2D([], [], marker="o", ls="", color="0.5", markersize=2, label="Modern (tiny dots)"))
    if len(mod_iran):
        h2.append(Line2D([], [], marker="s", ls="", color=col_of[mod_iran["sp"].iloc[0]], markeredgecolor="black",
                         markersize=float(np.sqrt(IRAN_SIZE)), label=f"Modern Iranians (n={len(mod_iran)})"))
    lax = fig.add_axes([0.73, 0.02, 0.26, 0.94]); lax.axis("off")
    l1 = lax.legend(handles=h1, loc="upper left", bbox_to_anchor=(0, 1), fontsize=8, frameon=False,
                    title="Colour: super-population", title_fontsize=9, labelspacing=0.45)
    l1._legend_box.align = "left"
    lax.add_artist(l1)
    l2 = lax.legend(handles=h2, loc="lower left", bbox_to_anchor=(0, 0), fontsize=8, frameon=False,
                    title="Size: age of ancient sample", title_fontsize=9, labelspacing=1.1)
    l2._legend_box.align = "left"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200)
    plt.close(fig)
    print("wrote", f"{out}.png", f"{out}.pdf")

if zoom:
    def tag(r, pre):
        return f"{pre}{r[0]:g}_{r[1]:g}" if r is not None else f"{pre}auto"
    plot_pair(pc_cols[0], pc_cols[1], f"{prefix}.pca_PC1_PC2.zoom_{tag(xr, 'x')}_{tag(yr, 'y')}", xr, yr)
else:
    plot_pair(pc_cols[0], pc_cols[1], f"{prefix}.pca_PC1_PC2")
    if len(pc_cols) >= 4:
        plot_pair(pc_cols[2], pc_cols[3], f"{prefix}.pca_PC3_PC4")

out = df[["sid", "sp", "is_ancient", "age"]].copy()
out["colour"] = out["sp"].map(col_of)
out.rename(columns={"sid": "IID", "sp": "super_population", "age": "age_BP"}).to_csv(
    f"{prefix}.pca_samples.tsv", sep="\t", index=False)
print("wrote", f"{prefix}.pca_samples.tsv")
PY
