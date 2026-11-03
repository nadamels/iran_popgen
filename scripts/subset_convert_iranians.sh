#!/usr/bin/env bash
# Subset modern Iranians from an EIGENSTRAT dataset (.geno/.ind/.snp) and write PLINK bed/bim/fam.
# Requires: convertf (EIGENSOFT) and plink (1.9) or plink2 on PATH.
# Usage: subset_convert_iranians.sh <eigenstrat_prefix> <out_prefix> [keep_list]
#   e.g. subset_convert_iranians.sh /path/to/v66.p1_compatibility_HO  modern_iran
set -euo pipefail
IN=${1:?eigenstrat prefix (without .geno/.ind/.snp)}
OUT=${2:?output prefix}
KEEP=${3:-"$(dirname "$0")/../metadata/modern_iranian_samples.txt"}

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

# 2) subset with plink (FID=IID=sample name, as familynames: NO)
awk '{print $1,$1}' "$KEEP" > "${OUT}.keep"
plink --bfile "${OUT}.full" --keep "${OUT}.keep" --make-bed --out "${OUT}"
rm -f "${OUT}.full".{bed,bim,fam} "${OUT}.convertf.par"
echo "Done: ${OUT}.bed/.bim/.fam"
