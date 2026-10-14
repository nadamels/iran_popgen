#!/bin/bash

#SBATCH -c 28
#SBATCH -t 6-00:00:00

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro

#Loading File as myenviro contains PLINK as well as most dependanceis for popgen.

# Rewrites the IDs of a PLINK binary fileset so that
#   new IID = current FID (column 1)   -> column 2 = old column 1
#   new FID = 0
# and writes the result to a new .bed/.bim/.fam set (so the .bed is regenerated
# consistently with the .fam). The original files are never modified.
#
# Usage: sbatch rename_fam_ids.sh [input_prefix] [output_prefix]
# Defaults apply to ad_anc_iran_merged.filtered (.bed/.bim/.fam)

set -euo pipefail

IN=${1:-ad_anc_iran_merged.filtered}
OUT=${2:-${IN}.renamed}
TMP=${IN}.update_ids.txt

for ext in bed bim fam; do
    [[ -s ${IN}.${ext} ]] || { echo "ERROR: ${IN}.${ext} not found" >&2; exit 1; }
done

# PLINK --update-ids format: old_FID old_IID new_FID new_IID
awk '{print $1, $2, 0, $1}' "${IN}.fam" > "${TMP}"

# Safety: new IIDs (old FIDs) must be unique, otherwise samples would collapse
# into duplicates (e.g. if every FID is already 0).
dups=$(cut -d' ' -f4 "${TMP}" | sort | uniq -d | wc -l)
if [[ ${dups} -gt 0 ]]; then
    echo "ERROR: ${dups} duplicated IID(s) would be created (old FIDs are not unique)." >&2
    echo "       Check column 1 of ${IN}.fam; nothing was written." >&2
    rm -f "${TMP}"
    exit 1
fi

plink --bfile "${IN}" \
      --update-ids "${TMP}" \
      --make-bed \
      --out "${OUT}"

# Verify: column 1 all 0, column 2 equals old column 1, rows in same order
paste -d' ' <(awk '{print $1}' "${IN}.fam") <(awk '{print $1, $2}' "${OUT}.fam") \
    | awk '$2!=0 || $3!=$1 {bad++} END{ if (bad) {print "ERROR: " bad " rows mismatched" > "/dev/stderr"; exit 1} else print "OK: all " NR " rows renamed" }'

rm -f "${TMP}"
