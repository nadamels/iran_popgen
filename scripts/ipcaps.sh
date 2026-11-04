#!/bin/bash

#SBATCH -c 10
#SBATCH -t 24:00:00
#SBATCH --mem=64G

# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/glimpse2

# IPCAPS (Iterative Pruning to Capture Population Structure; Chaichoompu et al. 2019, Source Code Biol
# Med 14:2, doi:10.1186/s13029-019-0072-6) on a PLINK fileset, restricted to the samples that are in BOTH
# the fileset and the mega sample table. IPCAPS is an unsupervised, iterative method: PCA, split the samples
# into sub-groups (rubikClust + Gaussian mixture modelling on the leading PCs), repeat inside every
# sub-group, and stop a branch when it has too few samples, the Fst between its sub-groups is too small,
# or EigenFit (difference of consecutive eigenvalues) falls below the threshold. Outliers are detected and
# removed at every iteration. The Population column of the table is used only to LABEL the results.
#   1. intersect fam IIDs with the table IIDs -> keep list + label file (same order as the subset fam)
#   2. subset the samples (plink --keep --make-bed): IPCAPS reads the whole fileset, so this is needed
#   3. ipcaps() in R
#   4. sample -> group table, Population-by-group table, per-group summary, optional top discriminators
#
# Usage:   sbatch ipcaps.sh <input_prefix> <mega_sample_table.xlsx|.csv> [output_dir]
#          bash   ipcaps.sh <input_prefix> <mega_sample_table.xlsx|.csv> [output_dir]
#   <input_prefix>   PLINK fileset without .bed/.bim/.fam. The paper's QC (see QC below) is assumed done;
#                    like umap.sh, the fileset is assumed to be LD-pruned already.
#   output_dir       default: <input_prefix>.ipcaps  (if it exists, IPCAPS writes to cluster_out, cluster_out1, ...
#                    inside it; the real location is printed at the end)
#
# Outputs (in the IPCAPS output directory):
#   groups.txt, tree_text.html, tree_scatter_cluster.html, tree_scatter_label.html, tree_scree.html, ...
#                                    IPCAPS's own results (group tree, PC plots coloured by cluster / label)
#   ipcaps_groups_with_IID.tsv       IID, group, node, label (Population) for every sample
#   ipcaps_label_by_group.tsv        Population x group counts
#   ipcaps_group_summary.tsv         per group: size, populations, top population and its fraction
#   top_discriminators_g<A>_g<B>.tsv top Fst SNPs between two groups (only when TOP_PAIRS is set)
#
# Settings (environment variables):
#   LABEL_COL="Population"   table column used to label the samples
#   METHOD=mix        internal clustering: mix (default) | mixmod | clara | pam | meanshift | apcluster | hclust
#   THRESHOLD=0.18    EigenFit cutoff, 0.03 - 0.18; smaller = potentially finer sub-structure (default 0.18)
#   MIN_FST=0.0008    minimum Fst between a pair of sub-groups
#   MIN_IN_GROUP=3    minimum number of individuals in a cluster; smaller clusters are counted as outliers
#                     (IPCAPS default is 20; 3 lets small populations form their own cluster)
#   TOP_PAIRS=        e.g. "4,5 4,6": run top.discriminator() (SNPs with the highest Fst) for these group pairs
#   NUM_TOP=100       SNPs reported per pair
#   AUTOSOME=0        1 = keep chromosomes 1-22 only when subsetting
#   QC=0              1 = also apply the paper's suggested PLINK QC when subsetting: --geno 0.02 --mind 0.05
#                     --maf 0.05 --hwe 0.001. Off by default: with ancient samples in the data, a 2% SNP
#                     missingness filter would remove most variants.
#   NO_PLOT=0         1 = no plots (for nodes without X11 / cairo graphics)
#
# Notes: IPCAPS loads the whole genotype matrix in memory (samples x SNPs) and missing genotypes are imputed
#   SNP-wise by the median. Run time grows quickly with the number of samples.
#
# Needs in the conda env: plink (1.9), python3 + pandas, R (>= 3.5) with the IPCAPS package, which also pulls
#   in KRIS, Rmixmod, fpc, LPCM, apcluster, expm, MASS and Matrix:
#     Rscript -e 'install.packages("IPCAPS")'    # or remotes::install_github("kridsadakorn/IPCAPS")
#   (openpyxl too if the table is .xlsx)

set -euo pipefail

if [ "${#ARGS[@]}" -lt 2 ] || [ "${#ARGS[@]}" -gt 3 ]; then
	echo "Usage: $0 <input_prefix> <mega_sample_table.xlsx|.csv> [output_dir]" >&2
	exit 1
fi

strip_ext() { local f=${1%.bed}; f=${f%.bim}; f=${f%.fam}; echo "$f"; }
abspath() { local d; d=$(cd "$(dirname "$1")" && pwd) && echo "${d}/$(basename "$1")"; }

IN=$(abspath "$(strip_ext "${ARGS[0]}")")
TABLE="${ARGS[1]}"
OUT="${ARGS[2]:-${IN}.ipcaps}"
mkdir -p "$(dirname "$OUT")"
OUT=$(abspath "$OUT")
LABEL_COL="${LABEL_COL:-Population}"
METHOD="${METHOD:-mix}"
THRESHOLD="${THRESHOLD:-0.18}"
MIN_FST="${MIN_FST:-0.0008}"
MIN_IN_GROUP="${MIN_IN_GROUP:-3}"
TOP_PAIRS="${TOP_PAIRS:-}"
NUM_TOP="${NUM_TOP:-100}"
AUTOSOME="${AUTOSOME:-0}"
QC="${QC:-0}"
NO_PLOT="${NO_PLOT:-0}"
THREADS=${SLURM_CPUS_PER_TASK:-${SLURM_JOB_CPUS_PER_NODE:-4}}

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
[[ "$METHOD" =~ ^(mix|mixmod|clara|pam|meanshift|apcluster|hclust)$ ]] || { echo "ERROR: unknown METHOD '$METHOD'" >&2; exit 1; }
[[ "$MIN_IN_GROUP" =~ ^[0-9]+$ ]] || { echo "ERROR: MIN_IN_GROUP must be an integer (got '$MIN_IN_GROUP')" >&2; exit 1; }
command -v plink >/dev/null || { echo "ERROR: plink not found in the active environment" >&2; exit 1; }
python3 -c "import pandas" 2>/dev/null || { echo "ERROR: python3 needs pandas" >&2; exit 1; }
command -v Rscript >/dev/null || { echo "ERROR: Rscript not found in the active environment (conda install -c conda-forge r-base)" >&2; exit 1; }
Rscript -e 'if (!requireNamespace("IPCAPS", quietly = TRUE)) quit(status = 1)' 2>/dev/null || {
	echo "ERROR: the R package IPCAPS is not installed: Rscript -e 'install.packages(\"IPCAPS\")'" >&2; exit 1; }

echo "input  : $IN  ($(wc -l < "${IN}.fam") samples, $(wc -l < "${IN}.bim") variants)"
echo "table  : $TABLE   label column: $LABEL_COL"
echo "ipcaps : method=$METHOD  threshold=$THRESHOLD  min.fst=$MIN_FST  min.in.group=$MIN_IN_GROUP   autosome=$AUTOSOME  qc=$QC"

TMP=$(mktemp -d "${OUT}.tmp.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# ---- 1. fam IIDs  ∩  table IIDs -> keep list, label file, ID list (all in fam order) -------------
python3 -I - "${IN}.fam" "$TABLE" "$LABEL_COL" "$TMP" <<'PY'
import re
import sys
import pandas as pd

fam_path, table_path, label_col, tmp = sys.argv[1:5]
fam = pd.read_csv(fam_path, sep=r"\s+", header=None, dtype=str, usecols=[0, 1], names=["FID", "IID"])
# Un-repaired fam where every IID is one constant string (ID stored in FID): fall back to FID.
use_fid = bool(fam["IID"].duplicated().any() and not fam["FID"].duplicated().any())
if use_fid:
    print("WARNING: IIDs in the fam are not unique but FIDs are; matching the table on FID")
fam["sid"] = (fam["FID"] if use_fid else fam["IID"]).str.strip()

if table_path.lower().endswith((".xlsx", ".xls")):
    tab = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
else:
    tab = pd.read_csv(table_path, dtype=str)
for c in ("IID", label_col):
    if c not in tab.columns:
        sys.exit(f"ERROR: table lacks column {c!r}; columns are: {', '.join(tab.columns)}")
tab["sid"] = tab["IID"].str.strip()
tab = tab.drop_duplicates("sid")

in_both = fam["sid"].isin(set(tab["sid"]))
print(f"fam samples: {len(fam)}   table samples: {len(tab)}   intersection: {int(in_both.sum())}")
if (~in_both).any():
    miss = fam.loc[~in_both, "sid"]
    print(f"WARNING: {len(miss)} fam sample(s) not in the table are dropped, e.g. {', '.join(miss.head(10))}")
print(f"{len(set(tab['sid']) - set(fam['sid']))} table sample(s) are not in this fileset (ignored)")
if in_both.sum() < 50:
    sys.exit("ERROR: fewer than 50 samples in common between the fam and the table")
kept = fam[in_both].reset_index(drop=True)             # fam order is preserved
if kept["sid"].duplicated().any():
    sys.exit("ERROR: duplicated IDs among the shared samples; run remove_duplicates.sh first")
kept = kept.merge(tab[["sid", label_col]], on="sid", how="left")

# IPCAPS reads the label file with read.table(header = FALSE): whitespace separated, so labels must not
# contain spaces, quotes or '#'.
def clean(x):
    x = "" if pd.isna(x) else str(x).strip()
    x = re.sub(r"[\s\"'#]+", "_", x)
    return x if x else "unlabelled"
kept["label"] = kept[label_col].map(clean)
kept[["FID", "IID"]].to_csv(f"{tmp}/keep.txt", sep="\t", header=False, index=False)
kept[["sid", "label"]].to_csv(f"{tmp}/labels.txt", sep=" ", header=False, index=False)
kept[["sid"]].to_csv(f"{tmp}/ids.txt", header=False, index=False)
print(f"{kept['label'].nunique()} distinct {label_col} labels among {len(kept)} samples")
PY

# ---- 2. subset the samples (IPCAPS reads the whole fileset) -----------------------------------
PLINK_OPTS=(--allow-no-sex --threads "$THREADS")
[ "$AUTOSOME" = "1" ] && PLINK_OPTS+=(--autosome)
[ "$QC" = "1" ] && PLINK_OPTS+=(--geno 0.02 --mind 0.05 --maf 0.05 --hwe 0.001)
plink --bfile "$IN" --keep "${TMP}/keep.txt" "${PLINK_OPTS[@]}" --make-bed --out "${TMP}/subset"
N_KEPT=$(wc -l < "${TMP}/keep.txt")
N_SUB=$(wc -l < "${TMP}/subset.fam")
echo "samples: $N_SUB of $N_KEPT shared samples kept   variants: $(wc -l < "${TMP}/subset.bim")"
if [ "$N_SUB" -ne "$N_KEPT" ]; then
	# --mind can drop samples; the label file is read by position, so rebuild it in the subset's fam order
	awk 'NR==FNR{lab[$1]=$2; next} {id = ($2 in lab) ? $2 : $1; print id, lab[id]}' "${TMP}/labels.txt" "${TMP}/subset.fam" > "${TMP}/labels2.txt"
	awk 'NR==FNR{lab[$1]=1; next} {print (($2 in lab) ? $2 : $1)}' "${TMP}/labels.txt" "${TMP}/subset.fam" > "${TMP}/ids.txt"
	mv "${TMP}/labels2.txt" "${TMP}/labels.txt"
fi

# ---- 3./4. IPCAPS and result tables (R) -------------------------------------------------------
cat > "${TMP}/run_ipcaps.R" <<'R'
args <- commandArgs(trailingOnly = TRUE)
bed <- paste0(args[1], ".bed"); label_file <- args[2]; ids_file <- args[3]; out_dir <- args[4]
method <- args[5]; threshold <- as.numeric(args[6]); min_fst <- as.numeric(args[7])
min_in_group <- as.integer(args[8]); no_plot <- args[9] == "1"; top_pairs <- args[10]; num_top <- as.integer(args[11])

suppressPackageStartupMessages(library(IPCAPS))
ids <- readLines(ids_file)
res <- ipcaps(bed = bed, label.file = label_file, lab.col = 2, out = out_dir, method = method,
              threshold = threshold, min.fst = min_fst, min.in.group = min_in_group, no.plot = no_plot)
if (is.null(res)) stop("ipcaps() returned NULL: the input files are not in the expected format")

res_dir <- res$output.dir
cl <- res$cluster                                   # columns: group, node, label, row.number
cl$IID <- ids[cl$row.number]
cl <- cl[order(cl$group, cl$IID), c("IID", "group", "node", "label")]
write.table(cl, file.path(res_dir, "ipcaps_groups_with_IID.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tab <- table(Population = cl$label, group = cl$group)
write.table(as.data.frame.matrix(tab), file.path(res_dir, "ipcaps_label_by_group.tsv"), sep = "\t", quote = FALSE,
            col.names = NA)

groups <- sort(unique(cl$group))
summ <- do.call(rbind, lapply(groups, function(g) {
  lab <- sort(table(cl$label[cl$group == g]), decreasing = TRUE)
  n <- sum(lab)
  data.frame(group = g, n = n, n_populations = length(lab), top_population = names(lab)[1],
             top_population_fraction = round(lab[[1]] / n, 3),
             populations = paste0(paste0(head(names(lab), 8), " (", head(as.integer(lab), 8), ")", collapse = "; "),
                                  if (length(lab) > 8) "; ..." else ""),
             small_group_outliers = n < min_in_group, stringsAsFactors = FALSE)
}))
write.table(summ, file.path(res_dir, "ipcaps_group_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n%d samples -> %d groups\n", nrow(cl), length(groups)))
print(summ[, c("group", "n", "n_populations", "top_population", "top_population_fraction")], row.names = FALSE)

pairs <- strsplit(trimws(top_pairs), "[[:space:]]+")[[1]]
for (p in pairs[nzchar(pairs)]) {
  g <- as.integer(strsplit(p, ",")[[1]])
  if (length(g) != 2 || anyNA(g)) { cat("NOTE: skipping malformed TOP_PAIRS entry", p, "\n"); next }
  top <- tryCatch(top.discriminator(res, g[1], g[2], num.top = num_top), error = function(e) {
    cat("NOTE: top.discriminator failed for groups", p, ":", conditionMessage(e), "\n"); NULL })
  if (!is.null(top))
    write.table(top, file.path(res_dir, sprintf("top_discriminators_g%d_g%d.tsv", g[1], g[2])),
                sep = "\t", quote = FALSE, row.names = FALSE)
}
cat("\nIPCAPS results are in:", res_dir, "\n")
R

Rscript "${TMP}/run_ipcaps.R" "${TMP}/subset" "${TMP}/labels.txt" "${TMP}/ids.txt" "$OUT" "$METHOD" "$THRESHOLD" \
	"$MIN_FST" "$MIN_IN_GROUP" "$NO_PLOT" "$TOP_PAIRS" "$NUM_TOP"
