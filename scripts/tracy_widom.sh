#!/bin/bash

#SBATCH -c 2
#SBATCH --mem=16G

# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/glimpse2

# Patterson-Price-Reich (2006) Tracy-Widom test for the number of informative principal components,
# the same test as EIGENSOFT smartpca / twstats (Patterson, Price & Reich 2006, PLoS Genet 2:e190).
# Restricted to the samples that are in BOTH the fileset and the mega sample table.
#   1. intersect fam IIDs with the table IIDs
#   2. read the genotypes straight from the .bed (no plink needed), normalise each SNP the smartpca way:
#      centre on the SNP mean, divide by sqrt(p(1-p)) with p = (1 + sum of genotypes) / (2 + 2 n), missing = 0
#   3. eigen-decomposition of the sample x sample covariance matrix (ALL eigenvalues)
#   4. for every eigenvalue i, use the eigenvalues i..end: rescale them to mean 1, estimate the effective
#      number of samples n' = m (m + 2) / (sum(l^2) - m), centre/scale the largest one with Johnstone's (2001)
#      formula and take the right tail of the Tracy-Widom (order 1) distribution -> p-value for PC i
#   5. informative PCs = the leading run of PCs with p < P_THRESHOLD (default 0.05)
# The input fileset is never modified. The SNPs should be LD-pruned and the samples unrelated (assumed).
#
# Usage:   sbatch tracy_widom.sh <input_prefix> <mega_sample_table.xlsx|.csv> [output_prefix]
#          bash   tracy_widom.sh <input_prefix> <mega_sample_table.xlsx|.csv> [output_prefix]
#   <input_prefix>   PLINK fileset without .bed/.bim/.fam (SNP-major .bed)
#   output_prefix    default: <input_prefix>.tw
#
# Outputs:  <out>.tw.tsv        PC, eigenvalue, % variance, effective n, TW statistic, p-value, informative
#           <out>.tw_pcs.tsv    IID, Population, PC1..PCk for the informative PCs (unit-length eigenvectors, the
#                               same scale as plink / smartpca eigenvec; at least 2 PCs are written)
#           <out>.tw.png/.pdf   scree plot and -log10(p) per PC
#
# Settings (environment variables):
#   P_THRESHOLD=0.05   a PC is significant when its p-value is below this
#   MIN_LEN=10         smartpca's minleneig: no test when fewer than this many eigenvalues remain
#   AUTOSOME=0         1 = use chromosomes 1-22 only
#
# Needs in the conda env: python3, numpy, scipy, pandas, matplotlib (openpyxl too if the table is .xlsx).
# Memory: about 8 bytes x SNPs x samples (361k SNPs x 68 samples is under 0.3 GB).

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
MIN_LEN="${MIN_LEN:-10}"
AUTOSOME="${AUTOSOME:-0}"

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
[[ "$MIN_LEN" =~ ^[0-9]+$ ]] || { echo "ERROR: MIN_LEN must be an integer (got '$MIN_LEN')" >&2; exit 1; }
python3 -c "import numpy, scipy, pandas, matplotlib" 2>/dev/null || {
	echo "ERROR: python3 needs numpy, scipy, pandas and matplotlib" >&2; exit 1; }

echo "input  : $IN  ($(wc -l < "${IN}.fam") samples, $(wc -l < "${IN}.bim") variants)"
echo "table  : $TABLE"
echo "output : ${OUT}.*   p-value threshold: $P_THRESHOLD   min eigenvalues left for a test: $MIN_LEN   autosomes only: $AUTOSOME"

python3 -I - "$IN" "$TABLE" "$OUT" "$P_THRESHOLD" "$MIN_LEN" "$AUTOSOME" <<'PY'
import os
import sys
import numpy as np
import pandas as pd
from scipy.special import airy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

prefix, table_path, out, pthr, min_len, autosome = sys.argv[1:7]
pthr, min_len, autosome = float(pthr), int(min_len), autosome == "1"


# ---- Tracy-Widom (order 1) right tail ------------------------------------------------------------
def tw1_tail(s):
    """P(TW1 > s) from the Fredholm determinant F1(s) = det(I - K) on L2(s, inf), K(x, y) = Ai((x + y) / 2) / 2
    (Bornemann 2010), by Gauss-Legendre quadrature. The tail is computed from the eigenvalues of K so that very
    small p-values keep their precision."""
    hi = max(s, 0.0) + 12.0
    n = int(min(300, 60 + 6 * (hi - s)))
    x, w = np.polynomial.legendre.leggauss(n)
    x = 0.5 * (hi - s) * x + 0.5 * (hi + s)
    w = 0.5 * (hi - s) * w
    K = 0.5 * airy(0.5 * (x[:, None] + x[None, :]))[0]
    mu = np.linalg.eigvalsh(np.sqrt(w)[:, None] * K * np.sqrt(w)[None, :])
    logF = np.log1p(-np.minimum(mu, 1 - 1e-15)).sum()
    return float(-np.expm1(logF))


def tw_norm(lam, p, n):
    """Johnstone (2001) centring and scaling, as twnorm() in EIGENSOFT."""
    if n < 0 or p < 0:
        return -10.0
    if n < p:
        return tw_norm(lam, n, p)
    y1 = np.sqrt(n - 1) + np.sqrt(p)
    return (lam - y1 * y1) / (y1 * (1 / np.sqrt(n - 1) + 1 / np.sqrt(p)) ** (1 / 3))


def tw_test(lam):
    """lam: eigenvalues in descending order, lam[0] is tested against the rest. EIGENSOFT dotwcalc()."""
    m = len(lam)
    ev = lam * (m / lam.sum())                       # rescale to mean 1
    bot = (ev ** 2).sum() - m
    if bot <= 0:
        return None
    zn = m * (m + 2) / bot                           # effective number of samples
    tw = tw_norm(ev[0] * zn, m, zn)
    return tw, tw1_tail(tw), zn


# ---- samples: fam intersect table -----------------------------------------------------------------
fam = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None, dtype=str, usecols=[0, 1], names=["FID", "IID"])
use_fid = bool(fam["IID"].duplicated().any() and not fam["FID"].duplicated().any())
if use_fid:
    print("WARNING: IIDs in the fam are not unique but FIDs are; matching the table on FID")
fam["sid"] = (fam["FID"] if use_fid else fam["IID"]).str.strip()
if table_path.lower().endswith((".xlsx", ".xls")):
    tab = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
else:
    tab = pd.read_csv(table_path, dtype=str)
if "IID" not in tab.columns:
    sys.exit(f"ERROR: {table_path} has no 'IID' column")
tab["sid"] = tab["IID"].str.strip()
tab = tab.drop_duplicates("sid")
sel = fam["sid"].isin(set(tab["sid"])).to_numpy()
print(f"fam samples: {len(fam)}   table samples: {len(tab)}   intersection: {int(sel.sum())}")
if (~sel).any():
    print(f"WARNING: {int((~sel).sum())} fam sample(s) not in the table are dropped, e.g. "
          + ", ".join(fam.loc[~sel, "sid"].head(10)))
n = int(sel.sum())
if n < 12:
    sys.exit(f"ERROR: only {n} shared samples; the test needs at least MIN_LEN + 2")
kept = fam[sel].reset_index(drop=True)
if kept["sid"].duplicated().any():
    sys.exit("ERROR: duplicated IDs among the shared samples; run remove_duplicates.sh first")

# ---- genotypes from the .bed ---------------------------------------------------------------------
bim = pd.read_csv(f"{prefix}.bim", sep=r"\s+", header=None, dtype=str, usecols=[0])
m_all, n_fam = len(bim), len(fam)
bps = (n_fam + 3) // 4
raw = np.fromfile(f"{prefix}.bed", dtype=np.uint8)
if raw.size < 3 or tuple(raw[:3]) != (0x6C, 0x1B, 0x01):
    sys.exit("ERROR: .bed is not a SNP-major PLINK 1 binary file (convert with plink --make-bed)")
if raw.size - 3 != m_all * bps:
    sys.exit(f"ERROR: .bed size does not match the .bim/.fam ({m_all} SNPs x {n_fam} samples)")
codes = ((raw[3:].reshape(m_all, bps)[:, :, None] >> np.array([0, 2, 4, 6], dtype=np.uint8)) & 3)
codes = codes.reshape(m_all, bps * 4)[:, :n_fam][:, sel]
del raw
G = np.array([2.0, np.nan, 1.0, 0.0], dtype=np.float32)[codes]          # 00 = hom A1, 01 = missing, 10 = het, 11 = hom A2
del codes
if autosome:
    chrom = pd.to_numeric(bim[0], errors="coerce").to_numpy()
    G = G[(chrom >= 1) & (chrom <= 22)]
nn = np.sum(~np.isnan(G), axis=1)
s = np.nansum(G, axis=1)
ok = nn >= 2
G, nn, s = G[ok], nn[ok], s[ok]
p = (1.0 + s) / (2.0 + 2.0 * nn)
X = (G - (s / nn)[:, None]) / np.sqrt(p * (1 - p))[:, None]
X = np.nan_to_num(X, nan=0.0).astype(np.float64)                          # missing genotypes -> 0 after centring
m = X.shape[0]
miss = 100 * np.mean(np.isnan(G))
del G
print(f"{n} samples, {m} SNPs used ({int((~ok).sum())} with fewer than 2 genotypes dropped); missing genotypes {miss:.2f}%")

# ---- eigen-decomposition and the test ------------------------------------------------------------
C = X.T @ X / m
del X
vals, vecs = np.linalg.eigh(C)
order = np.argsort(vals)[::-1]
vals, vecs = vals[order], vecs[:, order]
npos = int((vals > 1e-9 * vals[0]).sum())                                  # n - 1 non-trivial eigenvalues
lam = vals[:npos]
print(f"{npos} non-zero eigenvalues; leading eigenvalues: " + ", ".join(f"{v:.3f}" for v in lam[:5]))

rows = []
for i in range(npos):
    res = tw_test(lam[i:]) if npos - i >= min_len else None
    rows.append({"PC": i + 1, "eigenvalue": lam[i], "pct_variance": 100 * lam[i] / lam.sum(),
                 "effective_n": res[2] if res else np.nan, "TW_statistic": res[0] if res else np.nan,
                 "p_value": res[1] if res else np.nan})
tw = pd.DataFrame(rows)
tw["significant"] = tw["p_value"] < pthr
lead = 0
for sig in tw["significant"]:
    if not sig:
        break
    lead += 1
tw["informative"] = tw["PC"] <= lead
tw.to_csv(f"{out}.tw.tsv", sep="\t", index=False, float_format="%.6g")
print(tw.head(max(10, lead + 3)).to_string(index=False, float_format=lambda v: f"{v:.4g}"))
print(f"\n{int(tw['significant'].sum())} PCs have p < {pthr}; the leading run of significant PCs is PC1-PC{lead}"
      if lead else f"\nno PC has p < {pthr} (PC1 is not significant)")

# ---- informative PCs and plot --------------------------------------------------------------------
k = max(lead, 2)
pcs = pd.DataFrame(vecs[:, :k], columns=[f"PC{j + 1}" for j in range(k)])
pcs.insert(0, "IID", kept["sid"])
pop = kept[["sid"]].merge(tab[["sid"] + (["Population"] if "Population" in tab.columns else [])], on="sid", how="left")
if "Population" in pop.columns:
    pcs.insert(1, "Population", pop["Population"].to_numpy())
pcs.to_csv(f"{out}.tw_pcs.tsv", sep="\t", index=False, float_format="%.6g")

show = min(npos, max(20, lead + 5))
fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.5))
c = ["#d62728" if j < lead else "#7f7f7f" for j in range(show)]
a.bar(tw["PC"][:show], tw["pct_variance"][:show], color=c)
a.set_xlabel("PC")
a.set_ylabel("% of variance")
a.set_title(f"{os.path.basename(out)}: scree ({n} samples, {m} SNPs)")
pv = tw["p_value"][:show].clip(lower=1e-300)
b.bar(tw["PC"][:show], -np.log10(pv), color=c)
b.axhline(-np.log10(pthr), color="k", ls="--", lw=1, label=f"p = {pthr}")
b.set_xlabel("PC")
b.set_ylabel("-log10 p (Tracy-Widom)")
b.set_title(f"informative PCs (red): {lead}" if lead else "no significant PC")
b.legend(frameon=False)
fig.tight_layout()
fig.savefig(f"{out}.tw.png", dpi=200)
fig.savefig(f"{out}.tw.pdf")
print(f"wrote {out}.tw.tsv, {out}.tw_pcs.tsv, {out}.tw.png, {out}.tw.pdf")
PY
