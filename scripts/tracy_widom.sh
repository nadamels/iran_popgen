#!/bin/bash

#SBATCH -c 1
#SBATCH --mem=4G

# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/glimpse2

# Number of informative principal components: EIGENSOFT's twstats, i.e. the Patterson-Price-Reich (2006)
# Tracy-Widom test of each PC, run on the eigenvalues of a PCA you already have (PLINK 1.9 or PLINK 2).
#   1. read <prefix>.eigenval / <prefix>.eigenvec and check that ALL eigenvalues are there
#   2. twstats -i <eigenvalues> -o <out>   -> Tracy-Widom statistic and p-value for every PC
#   3. informative PCs = the leading run of PCs with p < P_THRESHOLD; table, PCs and plot
# twstats assumes independent markers: the PCA should come from LD-pruned SNPs (assumed).
#
# IMPORTANT: twstats needs ALL the eigenvalues. The test for PC i uses eigenvalues i, i+1, ..., last (the rest
# are the null "bulk" it is compared with), so a truncated list such as PLINK's default top 20 gives wrong
# p-values. With n samples run PLINK with  --pca <n-1>  (e.g. --pca 67 for 68 samples; the number of non-zero
# eigenvalues is n-1). The script stops if the eigenval file is shorter than that.
#
# Usage:   sbatch tracy_widom.sh <prefix> [mega_sample_table.xlsx|.csv] [output_prefix]
#          bash   tracy_widom.sh <prefix> [mega_sample_table.xlsx|.csv] [output_prefix]
#   <prefix>         expects  <prefix>.eigenvec  <prefix>.eigenval  (PLINK 1.9 or 2; eigenvec with or without a
#                    header line). <prefix>.bim, if present, is only used to count the SNPs
#   table            optional, only adds the Population column to the PCs file. Give "" to skip it
#   output_prefix    default: <prefix>.tw
#
# Outputs:  <out>.twstats.txt   twstats output as written by twstats
#           <out>.tw.tsv        PC, eigenvalue, % variance, TW statistic, p-value, effective n, informative
#           <out>.tw_pcs.tsv    IID, [Population,] PC1..PCk for the informative PCs (at least 2 PCs)
#           <out>.tw.png/.pdf   scree plot and -log10(p) per PC
#
# Settings (environment variables):
#   P_THRESHOLD=0.05   a PC is significant when its p-value is below this
#   MIN_LEN=10         twstats -m: NA when fewer than this many eigenvalues remain
#   TWTABLE=           optional path to EIGENSOFT's POPGEN/twtable (twstats -t); current twstats has the table built in
#
# Needs in the conda env: twstats from EIGENSOFT (conda install -c bioconda eigensoft, or
# https://github.com/DReichLab/EIG), python3 with numpy, pandas, matplotlib (openpyxl too if the table is .xlsx).

set -euo pipefail

if [ "${#ARGS[@]}" -lt 1 ] || [ "${#ARGS[@]}" -gt 3 ]; then
	echo "Usage: $0 <prefix> [mega_sample_table.xlsx|.csv] [output_prefix]" >&2
	exit 1
fi

abspath() { local d; d=$(cd "$(dirname "$1")" && pwd) && echo "${d}/$(basename "$1")"; }

PREFIX=$(abspath "${ARGS[0]%.eigenvec}")
PREFIX="${PREFIX%.eigenval}"
TABLE="${ARGS[1]:-}"
OUT="${ARGS[2]:-${PREFIX}.tw}"
mkdir -p "$(dirname "$OUT")"
OUT=$(abspath "$OUT")
P_THRESHOLD="${P_THRESHOLD:-0.05}"
MIN_LEN="${MIN_LEN:-10}"
TWTABLE="${TWTABLE:-}"

if [[ "$TABLE" == *.xlsx || "$TABLE" == *.xls ]] && ! python3 -c "import openpyxl" 2>/dev/null; then
	alt="${TABLE%.*}.csv"
	if [[ -s "$alt" ]]; then
		echo "NOTE: openpyxl not available; using $alt instead of $TABLE" >&2
		TABLE="$alt"
	fi
fi

for f in "${PREFIX}.eigenvec" "${PREFIX}.eigenval"; do
	[ -s "$f" ] || { echo "ERROR: $f not found or empty" >&2; exit 1; }
done
[ -z "$TABLE" ] || [ -s "$TABLE" ] || { echo "ERROR: table $TABLE not found" >&2; exit 1; }
[[ "$MIN_LEN" =~ ^[0-9]+$ ]] || { echo "ERROR: MIN_LEN must be an integer (got '$MIN_LEN')" >&2; exit 1; }
[ -z "$TWTABLE" ] || [ -s "$TWTABLE" ] || { echo "ERROR: TWTABLE $TWTABLE not found" >&2; exit 1; }
command -v twstats >/dev/null || { echo "ERROR: twstats not found in the active environment (conda install -c bioconda eigensoft)" >&2; exit 1; }
python3 -c "import numpy, pandas, matplotlib" 2>/dev/null || {
	echo "ERROR: python3 needs numpy, pandas and matplotlib" >&2; exit 1; }

# the eigenvec has a header line in PLINK 2 ("#FID IID PC1 ..." / "#IID PC1 ...") and none in PLINK 1.9
N_SAMPLES=$(grep -vc '^#' "${PREFIX}.eigenvec" || true)
N_EIG=$(grep -c '[0-9]' "${PREFIX}.eigenval" || true)
N_SNP=""
[ ! -s "${PREFIX}.bim" ] || N_SNP=$(wc -l < "${PREFIX}.bim")
echo "prefix : $PREFIX  ($N_SAMPLES samples, $N_EIG eigenvalues${N_SNP:+, $N_SNP variants})"
echo "output : ${OUT}.*   p-value threshold: $P_THRESHOLD"

# ---- 1. all eigenvalues present? ---------------------------------------------------------------
if [ "$N_EIG" -lt $((N_SAMPLES - 1)) ]; then
	echo "ERROR: ${PREFIX}.eigenval has $N_EIG eigenvalues but $N_SAMPLES samples have $((N_SAMPLES - 1)) non-zero ones." >&2
	echo "       twstats needs all of them (the p-value of PC i is computed against eigenvalues i..last); with only" >&2
	echo "       the top $N_EIG the p-values are wrong. Rerun the PCA with all of them, e.g. plink --pca $((N_SAMPLES - 1))" >&2
	exit 1
fi

# ---- 2. twstats --------------------------------------------------------------------------------
TW_ARGS=(-i "${PREFIX}.eigenval" -o "${OUT}.twstats.txt" -m "$MIN_LEN")
[ -z "$TWTABLE" ] || TW_ARGS+=(-t "$TWTABLE")
# twstats README: -n (#samples - 1) is needed when there are more samples than markers
if [ -n "$N_SNP" ] && [ "$N_SAMPLES" -gt "$N_SNP" ]; then TW_ARGS+=(-n $((N_SAMPLES - 1))); fi
twstats "${TW_ARGS[@]}"
[ -s "${OUT}.twstats.txt" ] || { echo "ERROR: twstats wrote no output" >&2; exit 1; }

# ---- 3. informative PCs, table, plot -----------------------------------------------------------
python3 -I - "${OUT}" "${PREFIX}.eigenvec" "$TABLE" "$P_THRESHOLD" "$N_SAMPLES" "${N_SNP:-NA}" <<'PY'
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out, evec_path, table_path, pthr, n_samples, n_snp = sys.argv[1:7]
pthr = float(pthr)

# twstats output: "#N  eigenvalue  difference  twstat  p-value  effect. n", NA where it is not computed
rows = []
for line in open(f"{out}.twstats.txt"):
    f = line.split()
    if not f or f[0].startswith("#"):
        continue
    num = [np.nan if v == "NA" else float(v) for v in f]
    rows.append({"PC": int(num[0]), "eigenvalue": num[1], "TW_statistic": num[3],
                 "p_value": num[4], "effective_n": num[5] if len(num) > 5 else np.nan})
if not rows:
    sys.exit(f"ERROR: could not read any rows from {out}.twstats.txt")
tw = pd.DataFrame(rows)
tw.insert(2, "pct_variance", 100 * tw["eigenvalue"] / tw["eigenvalue"].sum())
tw["significant"] = tw["p_value"] < pthr
lead = 0
for sig in tw["significant"]:
    if not sig:
        break
    lead += 1
tw["informative"] = tw["PC"] <= lead
tw.to_csv(f"{out}.tw.tsv", sep="\t", index=False, float_format="%.6g")
print(tw.head(max(10, lead + 3)).to_string(index=False, float_format=lambda v: f"{v:.4g}"))
if lead:
    print(f"\n{int(tw['significant'].sum())} PCs have p < {pthr}; the leading run of significant PCs is PC1-PC{lead}")
else:
    print(f"\nno PC has p < {pthr} (PC1 is not significant)")

# PLINK eigenvec: "#FID IID PC1 ..." / "#IID PC1 ..." (PLINK 2) or "FID IID PC1 ..." without header (PLINK 1.9)
with open(evec_path) as fh:
    first = fh.readline().split()
if first and first[0].startswith("#"):
    ev = pd.read_csv(evec_path, sep=r"\s+", dtype=str)
    ev.columns = [c.lstrip("#") for c in ev.columns]
    pc_cols = [c for c in ev.columns if c.upper().startswith("PC")]
else:
    ev = pd.read_csv(evec_path, sep=r"\s+", header=None, dtype=str)
    pc_cols = [f"PC{j}" for j in range(1, ev.shape[1] - 1)]
    ev.columns = ["FID", "IID"] + pc_cols
k = max(lead, 2)
if k > len(pc_cols):
    print(f"WARNING: the eigenvec has only {len(pc_cols)} PCs but {k} are informative; writing {len(pc_cols)}")
    k = len(pc_cols)
pcs = ev[["IID"] + pc_cols[:k]].copy()
pcs[pc_cols[:k]] = pcs[pc_cols[:k]].astype(float)
if table_path:
    if table_path.lower().endswith((".xlsx", ".xls")):
        tab = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
    else:
        tab = pd.read_csv(table_path, dtype=str)
    if "Population" in tab.columns and "IID" in tab.columns:
        tab["IID"] = tab["IID"].str.strip()
        pcs = pcs.merge(tab.drop_duplicates("IID")[["IID", "Population"]], on="IID", how="left")
        pcs = pcs[["IID", "Population"] + pc_cols[:k]]
        print(f"{int(pcs['Population'].notna().sum())} of {len(pcs)} samples have a Population in the table")
pcs.to_csv(f"{out}.tw_pcs.tsv", sep="\t", index=False, float_format="%.6g")

show = min(len(tw), max(20, lead + 5))
c = ["#d62728" if j < lead else "#7f7f7f" for j in range(show)]
fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.5))
a.bar(tw["PC"][:show], tw["pct_variance"][:show], color=c)
a.set_xlabel("PC")
a.set_ylabel("% of variance")
a.set_title(f"{os.path.basename(out)}: scree ({n_samples} samples, {n_snp} SNPs)")
b.bar(tw["PC"][:show], -np.log10(tw["p_value"][:show].clip(lower=1e-300)), color=c)
b.axhline(-np.log10(pthr), color="k", ls="--", lw=1, label=f"p = {pthr}")
b.set_xlabel("PC")
b.set_ylabel("-log10 p (twstats)")
b.set_title(f"informative PCs (red): {lead}" if lead else "no significant PC")
b.legend(frameon=False)
fig.tight_layout()
fig.savefig(f"{out}.tw.png", dpi=200)
fig.savefig(f"{out}.tw.pdf")
print(f"wrote {out}.twstats.txt, {out}.tw.tsv, {out}.tw_pcs.tsv, {out}.tw.png, {out}.tw.pdf")
PY
