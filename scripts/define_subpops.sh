#!/bin/bash

#SBATCH -c 10
#SBATCH -t 12:00:00
#SBATCH --mem=32G

# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/glimpse2

# Define sub-populations within the modern Iranian samples.
#   1. subset the keep list, QC (geno/mind/maf, autosomes), remove close relatives (KING), LD-prune
#   2. PCA, and ADMIXTURE K=2..KMAX (REPS random seeds per K, best run kept, 10-fold cross-validation)
#   3. subpop_consensus.py: GMM (k by BIC), HDBSCAN and ADMIXTURE hard calls on the same samples ->
#      final sub-populations (GMM), per-sample confidence, agreement between the methods
#   4. validate: pairwise Fst between the final sub-populations (plink --fst, Weir & Cockerham)
#
# NOTE: this pipeline was NOT derived from doi:10.1111/mec.13822 (the paper could not be read from the
# environment it was written in). Compare the QC/clustering choices below with the paper's methods and adjust.
#
# Usage:   sbatch define_subpops.sh <input_prefix> [keep_list] [output_dir] [mega_sample_table]
#          bash   define_subpops.sh <input_prefix> [keep_list] [output_dir] [mega_sample_table]
#   <input_prefix>  PLINK fileset (no .bed/.bim/.fam) that contains the modern Iranians
#   keep_list       FID IID or IID-only list; default metadata/modern_iranian_keep.txt
#   output_dir      default: subpops_out
#
# Settings (environment variables):
#   GENO=0.05  MIND=0.1  MAF=0.01     QC thresholds (relatively mild: ~100 samples, ~600k SNP array data)
#   KING=0.0884                       KING kinship cutoff (0.0884 = 3rd degree); 0 = skip. Needs plink2.
#   LD="200 25 0.4"                   --indep-pairwise window step r2 (plink 1.9 variant-count window = "50 5 0.2")
#   NPCS=10                           PCs used for clustering
#   KMAX=10  REPS=5  CV=10            ADMIXTURE settings (KMAX=0 skips ADMIXTURE)
#   MIN_SIZE=5  MIN_POST=0.9  ADMIX_MIN=0.6   see subpop_consensus.py
#
# Needs: plink (1.9), plink2 (only for KING), admixture, python3 with numpy, pandas, scikit-learn >= 1.3, matplotlib.

set -euo pipefail
set -- "${ARGS[@]}"

IN="${1:?input prefix}"
BASE="${SLURM_SUBMIT_DIR:-$PWD}"
KEEP="${2:-$BASE/metadata/modern_iranian_keep.txt}"
OUTDIR="${3:-subpops_out}"
TABLE="${4:-}"
GENO="${GENO:-0.05}"; MIND="${MIND:-0.1}"; MAF="${MAF:-0.01}"; KING="${KING:-0.0884}"
LD="${LD:-200 25 0.4}"; NPCS="${NPCS:-10}"
KMAX="${KMAX:-10}"; REPS="${REPS:-5}"; CV="${CV:-10}"
MIN_SIZE="${MIN_SIZE:-5}"; MIN_POST="${MIN_POST:-0.9}"; ADMIX_MIN="${ADMIX_MIN:-0.6}"
T="${SLURM_CPUS_PER_TASK:-4}"

for f in "$IN.bed" "$IN.bim" "$IN.fam" "$KEEP"; do [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done
mkdir -p "$OUTDIR"; O="$OUTDIR/modern_iran"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# keep list -> FID IID by looking the IIDs up in the fam (works for IID-only lists)
awk 'NR==FNR{k[$NF]=1; next} ($2 in k){print $1,$2}' "$KEEP" "$IN.fam" > "$O.keep"
echo "keep list: $(wc -l < "$KEEP") requested, $(wc -l < "$O.keep") found in $IN.fam"
[[ -s "$O.keep" ]] || { echo "ERROR: no overlap between keep list and fam" >&2; exit 1; }

# 1. QC
plink --bfile "$IN" --keep "$O.keep" --chr 1-22 --geno "$GENO" --mind "$MIND" --maf "$MAF" \
      --make-bed --out "$O.qc" --threads "$T"

# relatives
cp "$O.qc.bed" "$O.unrel.bed"; cp "$O.qc.bim" "$O.unrel.bim"; cp "$O.qc.fam" "$O.unrel.fam"
if [[ "$KING" != "0" ]]; then
    if command -v plink2 >/dev/null; then
        plink2 --bfile "$O.qc" --king-cutoff "$KING" --out "$O.king" --threads "$T"
        plink --bfile "$O.qc" --remove "$O.king.king.cutoff.out.id" --make-bed --out "$O.unrel" --threads "$T"
    else
        echo "WARNING: plink2 not found, relatives NOT removed (set KING=0 to silence)" >&2
    fi
fi

# LD pruning
plink --bfile "$O.unrel" --indep-pairwise $LD --out "$O.ld" --threads "$T"
plink --bfile "$O.unrel" --extract "$O.ld.prune.in" --make-bed --out "$O.pruned" --threads "$T"
echo "pruned dataset: $(wc -l < "$O.pruned.fam") samples, $(wc -l < "$O.pruned.bim") SNPs"

# 2. PCA
plink --bfile "$O.pruned" --pca 20 --out "$O.pca" --threads "$T"

# ADMIXTURE: REPS seeds per K, keep the run with the best log-likelihood; CV error from the same runs
QDIR="-"; CVT="-"
if [[ "$KMAX" -ge 2 ]]; then
    QDIR="$OUTDIR/admixture"; mkdir -p "$QDIR"
    PB="$(cd "$(dirname "$O.pruned")" && pwd)/$(basename "$O.pruned")"
    printf "K\tCV\n" > "$QDIR/cv_table.tsv"
    for K in $(seq 2 "$KMAX"); do
        best_ll=""; best=1; bestcv=""
        for r in $(seq 1 "$REPS"); do
            d="$QDIR/K${K}_r${r}"; mkdir -p "$d"
            (cd "$d" && admixture --cv="$CV" -s "$r" -j"$T" "$PB.bed" "$K" > log.txt 2>&1)
            ll=$(awk '/^Loglikelihood:/{print $2}' "$d/log.txt" | tail -1)
            cv=$(awk '/^CV error/{print $NF}' "$d/log.txt")
            if [[ -z "$best_ll" ]] || awk -v a="$ll" -v b="$best_ll" 'BEGIN{exit !(a>b)}'; then best_ll="$ll"; best=$r; bestcv="$cv"; fi
        done
        cp "$QDIR/K${K}_r${best}/$(basename "$PB").$K.Q" "$QDIR/modern_iran.$K.Q"
        printf "%s\t%s\n" "$K" "$bestcv" >> "$QDIR/cv_table.tsv"
        echo "K=$K best seed=$best loglik=$best_ll CV=$bestcv"
    done
    CVT="$QDIR/cv_table.tsv"
fi

# 3. consensus sub-populations
python3 -I "$HERE/subpop_consensus.py" "$O.pca.eigenvec" - "$CVT" "$QDIR" modern_iran "$O.pruned.fam" "$O" ${TABLE:+"$TABLE"} \
    --npcs "$NPCS" --min-size "$MIN_SIZE" --min-post "$MIN_POST" --admix-min "$ADMIX_MIN"

# 4. Fst between final sub-populations (pairs need >= 2 samples each; "unassigned" excluded)
awk 'NR==FNR{if(FNR>1 && $2!="unassigned") s[$1]=$2; next} ($2 in s){print $1,$2,s[$2]}' \
    "$O.subpops.tsv" "$O.pruned.fam" > "$O.within"
mapfile -t SPS < <(cut -d' ' -f3 "$O.within" | sort -u)
printf "pop1\tpop2\tmean_fst\tweighted_fst\n" > "$O.pairwise_fst.tsv"
for ((i=0; i<${#SPS[@]}; i++)); do
    for ((j=i+1; j<${#SPS[@]}; j++)); do
        a=${SPS[i]}; b=${SPS[j]}
        plink --bfile "$O.pruned" --within "$O.within" --keep-cluster-names "$a" "$b" --fst \
              --out "$O.fst_${a}_${b}" --threads "$T" > /dev/null || true
        m=$(awk '/Mean Fst estimate/{print $NF}' "$O.fst_${a}_${b}.log")
        w=$(awk '/Weighted Fst estimate/{print $NF}' "$O.fst_${a}_${b}.log")
        printf "%s\t%s\t%s\t%s\n" "$a" "$b" "$m" "$w" >> "$O.pairwise_fst.tsv"
    done
done
cat "$O.pairwise_fst.tsv"
echo "Done. Results: $O.subpops.tsv  $O.subpop_summary.tsv  $O.method_agreement.txt  $O.subpops_PCA.png"
