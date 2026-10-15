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

# Plot a PCA from the master sample table: COLOUR = country, SHAPE = period.
#
# Usage:   sbatch plot_pca.sh <prefix> [master_sample_table.tsv|.xlsx]
#          bash   plot_pca.sh <prefix> [master_sample_table.tsv|.xlsx]
#
#   <prefix>  expects  <prefix>.fam  <prefix>.eigenvec  <prefix>.eigenval
#             (PLINK 1.9 or PLINK 2 eigenvec; with or without a header line)
#   table     .tsv (preferred) or .xlsx (needs openpyxl; reads sheet "sample_table"); defaults to $MASTER_TABLE, else metadata/master_sample_table.tsv under
#             $SLURM_SUBMIT_DIR (or the current directory)
#
# Colour = country (colour per country is fixed by its rank among ALL countries in the master
#          table, so a country looks the same in every PCA from the same table).
# Shape  = archaeological period grouped into Mesolithic/Epipaleolithic, Neolithic, Chalcolithic,
#          Bronze Age, Iron Age, Historical (Anc/Rom/Byz/Medieval) - from archaeological_period,
#          or from date_mean_calBP when the period is missing.  Modern = small circle.
# Modern samples have no country in the table except the Iranians, so the rest are light grey.
# # Outputs (next to the prefix):  <prefix>.pca_PC1_PC2.{png,pdf}, <prefix>.pca_PC3_PC4.{png,pdf},
#                                <prefix>.pca_labels.tsv (country, colour, period, marker, n)

set -euo pipefail

PREFIX="${ARGS[0]:-}"
if [[ -z "$PREFIX" ]]; then
    echo "Usage: $0 <prefix> [master_sample_table.tsv|.xlsx]" >&2
    exit 1
fi
TABLE="${ARGS[1]:-${MASTER_TABLE:-${SLURM_SUBMIT_DIR:-$PWD}/metadata/master_sample_table.tsv}}"

for f in "${PREFIX}.fam" "${PREFIX}.eigenvec" "${PREFIX}.eigenval" "$TABLE"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done

python3 -c "import numpy, pandas, matplotlib" 2>/dev/null || {
    echo "ERROR: python3 in the active environment needs numpy, pandas and matplotlib" >&2; exit 1; }

echo "prefix : $PREFIX"
echo "table  : $TABLE"

python3 - "$PREFIX" "$TABLE" <<'PY'
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

prefix, table_path = sys.argv[1], sys.argv[2]

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

# ---- country + period from the master table ---------------------------------------------------
if table_path.lower().endswith((".xlsx", ".xls")):
    try:
        t = pd.read_excel(table_path, sheet_name="sample_table", dtype=str)
    except ImportError:
        sys.exit("ERROR: reading .xlsx needs openpyxl; pass master_sample_table.tsv instead")
else:
    t = pd.read_csv(table_path, sep="\t", dtype=str)
t = t.drop_duplicates("sample_id")

# Period -> broad archaeological class (marker shape).  Order = legend order.
PERIODS = ["Mesolithic/Epipaleolithic", "Neolithic", "Chalcolithic", "Bronze Age", "Iron Age",
           "Historical", "Unknown period"]
PERIOD_MAP = {
    "natufian": "Mesolithic/Epipaleolithic", "mes": "Mesolithic/Epipaleolithic",
    "ppn": "Neolithic", "n": "Neolithic", "en": "Neolithic", "m / n": "Neolithic",
    "chl": "Chalcolithic", "c / en": "Chalcolithic",
    "ba": "Bronze Age", "emba": "Bronze Age", "mlba": "Bronze Age", "lba": "Bronze Age",
    "ia": "Iron Age", "lba / ia": "Iron Age", "urartian": "Iron Age",
    "anc": "Historical", "rom": "Historical", "byz": "Historical", "romByz".lower(): "Historical",
    "mdv": "Historical", "h": "Historical",
}
def period_from_date(bp):                      # fallback only when the table gives no period
    if not np.isfinite(bp):
        return "Unknown period"
    for lim, name in ((2200, "Historical"), (3200, "Iron Age"), (5300, "Bronze Age"),
                      (7000, "Chalcolithic"), (11700, "Neolithic")):
        if bp < lim:
            return name
    return "Mesolithic/Epipaleolithic"

bp = pd.to_numeric(t["date_mean_calBP"], errors="coerce")
per = t["archaeological_period"].fillna("").str.strip().str.lower().map(PERIOD_MAP)
from_date = per.isna() & (t["sample_type"] == "ancient")
per = per.where(~from_date, bp.map(period_from_date))
t["period"] = np.where(t["sample_type"] == "ancient", per.fillna("Unknown period"), "Modern")
t["country_plot"] = t["country"].fillna("").str.strip().replace("", np.nan)

df = ev.merge(t[["sample_id", "sample_type", "country_plot", "period"]],
              left_on="sid", right_on="sample_id", how="left")
unk = df["sample_type"].isna()
if unk.any():
    print(f"WARNING: {unk.sum()} samples not in the master table (plotted grey 'x'): "
          f"{', '.join(df.loc[unk, 'sid'].head(10))}{' ...' if unk.sum() > 10 else ''}")
df.loc[unk, "period"] = "Not in table"
anc = df["sample_type"] == "ancient"
n_nc = int((anc & df["country_plot"].isna()).sum())
if n_nc:
    print(f"NOTE: {n_nc} ancient samples have no country; plotted as 'Country unknown'")
n_pd = int((anc & df["sample_id"].isin(t.loc[from_date, "sample_id"])).sum())
if n_pd:
    print(f"NOTE: {n_pd} ancient samples have no archaeological_period; shape inferred from date_mean_calBP")
df["country_plot"] = df["country_plot"].fillna("Country unknown")

# ---- stable colour per country, marker per period ----------------------------------------------
# Rank over ALL countries in the table so a country keeps its colour across plots.
countries = sorted(set(t["country_plot"].dropna()))
cols = np.vstack([plt.cm.tab20(np.arange(20)), plt.cm.tab20b(np.arange(20))])
col_of = {c: cols[i % len(cols)] for i, c in enumerate(countries)}
col_of["Country unknown"] = np.array([0.35, 0.35, 0.35, 1])
MODERN_GREY = np.array([0.82, 0.82, 0.82, 1])
MARK = dict(zip(PERIODS, ["P", "s", "^", "D", "v", "X", "*"]))
MARK["Modern"] = "o"

# ---- plot --------------------------------------------------------------------------------------
rows = []
def plot_pair(a, b, out):
    fig = plt.figure(figsize=(17, 10))
    ax = fig.add_axes([0.06, 0.07, 0.64, 0.88])
    # modern without country: grey backdrop; modern with country: small coloured circles
    m = df[df["period"] == "Modern"]
    g = m[m["country_plot"] == "Country unknown"]
    ax.scatter(g[a], g[b], s=12, c=[MODERN_GREY], marker="o", edgecolors="none", alpha=0.7, zorder=1)
    for ctry, s in m[m["country_plot"] != "Country unknown"].groupby("country_plot"):
        ax.scatter(s[a], s[b], s=18, c=[col_of[ctry]], marker="o", edgecolors="none", alpha=0.8, zorder=2)
    for per_name in PERIODS:
        for ctry, s in df[(df["sample_type"] == "ancient") & (df["period"] == per_name)].groupby("country_plot"):
            ax.scatter(s[a], s[b], s=48 if per_name != "Mesolithic/Epipaleolithic" else 60,
                       c=[col_of[ctry]], marker=MARK[per_name], edgecolors="black", linewidths=0.5,
                       alpha=0.95, zorder=4)
    u = df[df["period"] == "Not in table"]
    if len(u):
        ax.scatter(u[a], u[b], s=20, c="grey", marker="x", zorder=3)
    ia, ib = pc_cols.index(a), pc_cols.index(b)
    ax.set_xlabel(f"{a} ({pct[ia]:.2f}% of listed eigenvalues)" if ia < len(pct) else a)
    ax.set_ylabel(f"{b} ({pct[ib]:.2f}% of listed eigenvalues)" if ib < len(pct) else b)
    ax.set_title(f"{prefix.split('/')[-1]}  |  n={len(df)}  "
                 f"(ancient {int(anc.sum())}, modern {int((df.period=='Modern').sum())})")
    ax.grid(alpha=0.2)

    # legend 1: colour = country (ancient + modern with a country)
    used = df[(df["period"] != "Modern") | (df["country_plot"] != "Country unknown")]
    ctr = sorted(used["country_plot"].unique(), key=lambda c: (c == "Country unknown", c))
    h1 = [Line2D([], [], marker="o", ls="", color=col_of[c], markersize=7, markeredgecolor="black",
                 markeredgewidth=0.4, label=f"{c} ({int((used.country_plot == c).sum())})") for c in ctr]
    n_mg = int(((df.period == "Modern") & (df.country_plot == "Country unknown")).sum())
    if n_mg:
        h1.append(Line2D([], [], marker="o", ls="", color=MODERN_GREY, markersize=7,
                         label=f"Modern, no country ({n_mg})"))
    # legend 2: shape = period
    h2 = []
    for per_name in PERIODS + ["Modern"]:
        n = int((df.period == per_name).sum())
        if n:
            h2.append(Line2D([], [], marker=MARK[per_name], ls="", color="0.4", markersize=8,
                             markeredgecolor="black", markeredgewidth=0.4,
                             label=f"{per_name if per_name != 'Modern' else 'Modern (present-day)'} ({n})"))
    lax = fig.add_axes([0.72, 0.02, 0.27, 0.94]); lax.axis("off")
    l1 = lax.legend(handles=h1, loc="upper left", bbox_to_anchor=(0, 1), ncol=1, fontsize=8,
                    frameon=False, title="Colour: country", title_fontsize=9, labelspacing=0.3)
    l1._legend_box.align = "left"
    lax.add_artist(l1)
    l2 = lax.legend(handles=h2, loc="lower left", bbox_to_anchor=(0, 0), fontsize=8, frameon=False,
                    title="Shape: period", title_fontsize=9, labelspacing=0.4)
    l2._legend_box.align = "left"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200)
    plt.close(fig)
    print("wrote", f"{out}.png", f"{out}.pdf")

plot_pair(pc_cols[0], pc_cols[1], f"{prefix}.pca_PC1_PC2")
if len(pc_cols) >= 4:
    plot_pair(pc_cols[2], pc_cols[3], f"{prefix}.pca_PC3_PC4")

# key for reproducibility
for (per_name, ctry), s in df.groupby(["period", "country_plot"]):
    c = MODERN_GREY if (per_name == "Modern" and ctry == "Country unknown") else col_of[ctry]
    rows.append((ctry, matplotlib.colors.to_hex(c), per_name, MARK.get(per_name, "x"), len(s)))
pd.DataFrame(rows, columns=["country", "colour", "period", "marker", "n"]).to_csv(
    f"{prefix}.pca_labels.tsv", sep="\t", index=False)
print("wrote", f"{prefix}.pca_labels.tsv")
PY
