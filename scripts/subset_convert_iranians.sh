#!/bin/bash
#SBATCH -c 2

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro

# Subset modern Iranians from an EIGENSTRAT dataset (.geno/.ind/.snp, one shared prefix)
# and write PLINK bed/bim/fam. Requires convertf (EIGENSOFT) and plink in myenviro.
#
# Usage: sbatch subset_convert_iranians.sh <eigenstrat_prefix> <sample_list> <out_dir> <out_prefix>
#   <eigenstrat_prefix> full path without .geno/.ind/.snp
#   <sample_list>       1-column file of sample IDs to keep (matching the .ind IDs)
#   <out_dir>           output directory (created if missing)
#   <out_prefix>        prefix for the output bed/bim/fam
set -euo pipefail

IN=${1:?eigenstrat prefix (no extension)}
KEEP=${2:?1-column sample list}
OUTDIR=${3:?output directory}
PREFIX=${4:?output prefix}

for e in geno ind snp; do [[ -s "${IN}.${e}" ]] || { echo "Missing ${IN}.${e}" >&2; exit 1; }; done
[[ -s "$KEEP" ]] || { echo "Missing sample list $KEEP" >&2; exit 1; }
mkdir -p "$OUTDIR"
OUT="${OUTDIR}/${PREFIX}"

# 1) convertf: EIGENSTRAT -> PACKEDPED (all samples)
cat > "${OUT}.convertf.par" <<PAR
genotypename:    ${IN}.geno
snpname:         ${IN}.snp
indivname:       ${IN}.ind
outputformat:    PACKEDPED
genotypeoutname: ${OUT}.full.bed
snpoutname:      ${OUT}.full.bim
indivoutname:    ${OUT}.full.fam
familynames:     NO
PAR
convertf -p "${OUT}.convertf.par"

# 2) subset with plink (FID=IID=sample name because familynames: NO)
awk '{print $1,$1}' "$KEEP" > "${OUT}.keep"
plink --bfile "${OUT}.full" --keep "${OUT}.keep" --make-bed --out "${OUT}"

n=$(wc -l < "${OUT}.fam"); echo "Kept $n of $(wc -l < "$KEEP") requested samples"
rm -f "${OUT}.full".{bed,bim,fam} "${OUT}.convertf.par" "${OUT}.keep"
echo "Done: ${OUT}.bed/.bim/.fam"
