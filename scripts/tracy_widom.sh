#!/bin/bash

#SBATCH -c 4
#SBATCH --mem=16G

# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/glimpse2

# Number of informative principal components with EIGENSOFT: smartpca for the PCA, then twstats for the
# Patterson-Price-Reich (2006) Tracy-Widom test of each PC. Restricted to the samples that are in BOTH the
# PLINK fileset and the mega sample table.
#   1. intersect fam IIDs with the table IIDs; the other samples are set to "Ignore" in a copy of the fam
#      (EIGENSOFT's .pedind format), so the fileset itself is read as it is (PACKEDPED: .bed + .bim)
#   2. smartpca -p <parfile>   -> eigenvectors and ALL eigenvalues (no outlier removal)
#   3. twstats -i <eigenvalues> -o <out>   -> Tracy-Widom statistic and p-value for every PC
#   4. informative PCs = the leading run of PCs with p < P_THRESHOLD; table, PCs and plot
# The input fileset is never modified. twstats assumes independent markers: use LD-pruned SNPs (assumed).
#
# Usage:   sbatch tracy_widom.sh <input_prefix> <mega_sample_table.xlsx|.csv> [output_prefix]
#          bash   tracy_widom.sh <input_prefix> <mega_sample_table.xlsx|.csv> [output_prefix]
#   <input_prefix>   PLINK fileset without .bed/.bim/.fam (SNP-major .bed, SNPs in genome order)
#   output_prefix    default: <input_prefix>.tw
#
# Outputs:  <out>.twstats.txt   twstats output as written by twstats
#           <out>.tw.tsv        PC, eigenvalue, % variance, TW statistic, p-value, effective n, informative
#           <out>.tw_pcs.tsv    IID, Population, PC1..PCk for the informative PCs (smartpca eigenvectors; at least 2 PCs)
#           <out>.tw.png/.pdf   scree plot and -log10(p) per PC
#           <out>.smartpca.log  smartpca log
#
# Settings (environment variables):
#   P_THRESHOLD=0.05   a PC is significant when its p-value is below this
#   NPCS=20            eigenvectors written by smartpca (the test itself uses all eigenvalues). If every one
#                      of them is significant, raise it
#   MIN_LEN=10         twstats -m: NA when fewer than this many eigenvalues remain
#   TWTABLE=           optional path to EIGENSOFT's POPGEN/twtable (twstats -t); current twstats has the table built in
#
# smartpca defaults kept: SNPs on chromosome X are excluded (noxdata YES) and SNPs are normalised by allele
# frequency; outlier removal is switched off (numoutlieriter 0) so that all samples are used.
#
# Needs in the conda env: smartpca and twstats from EIGENSOFT (conda install -c bioconda eigensoft, or
# https://github.com/DReichLab/EIG), python3 with numpy, pandas, matplotlib (openpyxl too if the table is .xlsx).

set -euo pipefail

if [ "${#ARGS[@]}" -lt 2 ] || [ "${#ARGS[@]}" -gt 3 ]; then
	echo "Usage: $0 <input_prefix> <mega_sample_table.xlsx|.csv> [output_prefix]" >&2
	exit 1
fi

strip_ext() { local f=${1%.bed}; f=${f%.bim}; f=${f%.fam}; echo "$f"; }
abspath() { local d; d=$(cd "$(dirname "$1")" && pwd) && echo "${d}/$(basename "$1")"; }

IN=$(abspath "$(strip_ext "${ARGS[0]}")")
TABLE="${ARGS[1]}"
OUT="${ARGS[2]:-${IN}.tw}"
mkdir -p "$(dirname "$OUT")"
OUT=$(abspath "$OUT")
P_THRESHOLD="${P_THRESHOLD:-0.05}"
NPCS="${NPCS:-20}"
MIN_LEN="${MIN_LEN:-10}"
TWTABLE="${TWTABLE:-}"

if [[ "$TABLE" == *.xlsx || "$TABLE" == *.xls ]] && ! python3 -c "import openpyxl" 2>/dev/null; then
	alt="${TABLE%.*}.csv"
	if [[ -s "$alt" ]]; then
		echo "NOTE: openpyxl not available; using $alt instead of $TABLE" >&2
		TABLE="$alt"
	fi
fi

for ext in bed bim fam; do
	[ -s "${IN}.${ext}" ] || { echo "ERROR: ${IN}.${ext} not found" >&2; exit 1; }
done
[ -s "$TABLE" ] || { echo "ERROR: table $TABLE not found" >&2; exit 1; }
[[ "$NPCS" =~ ^[0-9]+$ && "$NPCS" -ge 2 ]] || { echo "ERROR: NPCS must be an integer >= 2 (got '$NPCS')" >&2; exit 1; }
[[ "$MIN_LEN" =~ ^[0-9]+$ ]] || { echo "ERROR: MIN_LEN must be an integer (got '$MIN_LEN')" >&2; exit 1; }
[ -z "$TWTABLE" ] || [ -s "$TWTABLE" ] || { echo "ERROR: TWTABLE $TWTABLE not found" >&2; exit 1; }
for tool in smartpca twstats; do
	command -v "$tool" >/dev/null || { echo "ERROR: $tool not found in the active environment (conda install -c bioconda eigensoft)" >&2; exit 1; }
done
python3 -c "import numpy, pandas, matplotlib" 2>/dev/null || {
	echo "ERROR: python3 needs numpy, pandas and matplotlib" >&2; exit 1; }

N_SNP=$(wc -l < "${IN}.bim")
echo "input  : $IN  ($(wc -l < "${IN}.fam") samples, $N_SNP variants)"
echo "table  : $TABLE"
echo "output : ${OUT}.*   p-value threshold: $P_THRESHOLD   eigenvectors written: $NPCS"

TMP=$(mktemp -d "${OUT}.tmp.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# ---- 1. fam IIDs  ∩  table IIDs -> indiv file with the other samples set to Ignore ---------------
python3 -I - "${IN}.fam" "$TABLE" "$TMP" <<'PY'
import sys
import pandas as pd

fam_path, table_path, tmp = sys.argv[1:4]
fam = pd.read_csv(fam_path, sep=r"\s+", header=None, dtype=str,
                  names=["FID", "IID", "PID", "MID", "SEX", "PHENO"])
fam["sid"] = fam["IID"].str.strip()
if table_path.lower().endswith((".xlsx", ".xls")):
    tab = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
else:
    tab = pd.read_csv(table_path, dtype=str)
if "IID" not in tab.columns:
    sys.exit(f"ERROR: {table_path} has no 'IID' column")
ids = set(tab["IID"].dropna().str.strip())

in_both = fam["sid"].isin(ids)
print(f"fam samples: {len(fam)}   table samples: {len(ids)}   intersection: {int(in_both.sum())}")
if (~in_both).any():
    print(f"WARNING: {int((~in_both).sum())} fam sample(s) not in the table are ignored, e.g. "
          + ", ".join(fam.loc[~in_both, "sid"].head(10)))
print(f"{len(ids - set(fam['sid']))} table sample(s) are not in this fileset (ignored)")
if in_both.sum() < 12:
    sys.exit("ERROR: fewer than 12 samples in common between the fam and the table")
if fam.loc[in_both, "sid"].duplicated().any():
    sys.exit("ERROR: duplicated IIDs among the shared samples (smartpca identifies samples by IID here); "
             "run remove_duplicates.sh first")
if (fam.loc[in_both, "sid"].str.len() > 39).any():
    sys.exit("ERROR: an IID is longer than 39 characters, the smartpca limit")

# .pedind: the fam columns; the 6th column is the label, "Ignore" removes the sample
ind = fam[["FID", "IID", "PID", "MID", "SEX"]].copy()
ind["LABEL"] = in_both.map({True: "Case", False: "Ignore"})
ind.to_csv(f"{tmp}/ind.pedind", sep=" ", header=False, index=False)
print(f"{int(in_both.sum())} samples are used")
PY
N_USED=$(awk '$6 == "Case"' "${TMP}/ind.pedind" | wc -l)

# ---- 2. smartpca (PACKEDPED: .bed + .bim + .pedind) ------------------------------------------
cat > "${TMP}/par" <<EOF
genotypename:    ${IN}.bed
snpname:         ${IN}.bim
indivname:       ${TMP}/ind.pedind
evecoutname:     ${TMP}/pca.evec
evaloutname:     ${TMP}/pca.eval
numoutevec:      ${NPCS}
numoutlieriter:  0
familynames:     NO
EOF
smartpca -p "${TMP}/par" > "${OUT}.smartpca.log" 2>&1 || {
	echo "ERROR: smartpca failed; see ${OUT}.smartpca.log" >&2; tail -n 15 "${OUT}.smartpca.log" >&2; exit 1; }
[ -s "${TMP}/pca.eval" ] || { echo "ERROR: smartpca wrote no eigenvalues; see ${OUT}.smartpca.log" >&2; exit 1; }
echo "smartpca: done ($(wc -l < "${TMP}/pca.eval") eigenvalues); log: ${OUT}.smartpca.log"

# ---- 3. twstats --------------------------------------------------------------------------------
TW_ARGS=(-i "${TMP}/pca.eval" -o "${OUT}.twstats.txt" -m "$MIN_LEN")
[ -z "$TWTABLE" ] || TW_ARGS+=(-t "$TWTABLE")
# twstats README: -n (#samples - 1) is needed when there are more samples than markers
[ "$N_USED" -le "$N_SNP" ] || TW_ARGS+=(-n $((N_USED - 1)))
twstats "${TW_ARGS[@]}"
[ -s "${OUT}.twstats.txt" ] || { echo "ERROR: twstats wrote no output" >&2; exit 1; }

# ---- 4. informative PCs, table, plot -----------------------------------------------------------
python3 -I - "${OUT}" "${TMP}/pca.evec" "$TABLE" "$P_THRESHOLD" "$NPCS" "$N_USED" "$N_SNP" <<'PY'
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out, evec_path, table_path, pthr, npcs, n_used, n_snp = sys.argv[1:8]
pthr, npcs = float(pthr), int(npcs)

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
if lead >= npcs:
    print(f"WARNING: all {npcs} PCs written by smartpca are significant; rerun with a larger NPCS")

# smartpca evec: "#eigvals: ..." header, then  ID  PC1 ... PCk  label
ev = pd.read_csv(evec_path, sep=r"\s+", header=None, skiprows=1, dtype=str)
k_ev = ev.shape[1] - 2
ev.columns = ["IID"] + [f"PC{j + 1}" for j in range(k_ev)] + ["label"]
k = min(max(lead, 2), k_ev)
pcs = ev[["IID"] + [f"PC{j + 1}" for j in range(k)]].copy()
if table_path.lower().endswith((".xlsx", ".xls")):
    tab = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
else:
    tab = pd.read_csv(table_path, dtype=str)
if "Population" in tab.columns:
    tab["IID"] = tab["IID"].str.strip()
    pcs = pcs.merge(tab.drop_duplicates("IID")[["IID", "Population"]], on="IID", how="left")
    pcs = pcs[["IID", "Population"] + [f"PC{j + 1}" for j in range(k)]]
pcs.to_csv(f"{out}.tw_pcs.tsv", sep="\t", index=False)

show = min(len(tw), max(20, lead + 5))
c = ["#d62728" if j < lead else "#7f7f7f" for j in range(show)]
fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.5))
a.bar(tw["PC"][:show], tw["pct_variance"][:show], color=c)
a.set_xlabel("PC")
a.set_ylabel("% of variance")
a.set_title(f"{os.path.basename(out)}: scree ({n_used} samples, {n_snp} SNPs)")
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
