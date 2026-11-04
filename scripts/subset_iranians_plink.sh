#!/bin/bash
#SBATCH -c 2

# `source`d scripts inherit our arguments and conda's activate rejects them (ArgumentError): hide, then restore
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro
set -- "${ARGS[@]}"

# Subsets a PLINK binary fileset (.bed/.bim/.fam) to the samples in a 1-column list.
# In the Reich-converted .fam, column 1 (FID) is just a row number and column 2 (IID)
# is the sample name, so --keep needs the real FID/IID pairs: they are looked up
# in the .fam by matching the list against column 2 (IID).
#
# Usage: sbatch subset_iranians_plink.sh <bfile_prefix> <sample_list> <out_dir> <out_prefix>
#   <bfile_prefix> full path without .bed/.bim/.fam
#   <sample_list>  1-column file of sample IDs (IID, e.g. metadata/modern_iranian_samples.txt)
#   <out_dir>      output directory (created if missing)
#   <out_prefix>   prefix for the output bed/bim/fam

set -euo pipefail

BFILE=${1:?bfile prefix (no extension)}
LIST=${2:?1-column sample list}
OUTDIR=${3:?output directory}
PREFIX=${4:?output prefix}

for e in bed bim fam; do [[ -s "${BFILE}.${e}" ]] || { echo "Missing ${BFILE}.${e}" >&2; exit 1; }; done
[[ -s "$LIST" ]] || { echo "Missing sample list $LIST" >&2; exit 1; }
mkdir -p "$OUTDIR"
OUT="${OUTDIR}/${PREFIX}"

# Build 2-column (FID IID) keep file by matching the list to .fam column 2
awk 'NR==FNR {want[$1]=1; next} ($2 in want) {print $1, $2}' "$LIST" "${BFILE}.fam" > "${OUT}.keep"

nwant=$(sort -u "$LIST" | wc -l)
nfound=$(wc -l < "${OUT}.keep")
echo "Requested $nwant samples; found $nfound in ${BFILE}.fam"
if [[ "$nfound" -lt "$nwant" ]]; then
    echo "WARNING: samples not found in .fam:" >&2
    awk 'NR==FNR {have[$2]=1; next} !($1 in have)' "${BFILE}.fam" "$LIST" >&2
fi
[[ "$nfound" -gt 0 ]] || { echo "No samples matched; aborting" >&2; exit 1; }

plink --bfile "$BFILE" --keep "${OUT}.keep" --make-bed --out "$OUT"

echo "Kept $(wc -l < "${OUT}.fam") samples -> ${OUT}.bed/.bim/.fam"
