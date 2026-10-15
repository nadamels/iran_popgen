#!/bin/bash

#SBATCH -c 10
#SBATCH -t 10:00:00
#SBATCH --mail-user=nme3923@nyu.edu
#SBATCH --mail-type=END,FAIL

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro

#Loading File
# Usage: sbatch merging.sh <file1_prefix> <file2_prefix> <output_prefix>
# Each prefix may be a relative or absolute path to a PLINK fileset
# (without .bed/.bim/.fam), and each can live in a different directory.
# The output prefix's directory is created if needed.

set -euo pipefail

if [ "$#" -ne 3 ]; then
	echo "Usage: $0 <file1_prefix> <file2_prefix> <output_prefix>" >&2
	exit 1
fi

# Resolve to absolute paths (strip an accidental .bed/.bim/.fam extension)
strip_ext() { local f=${1%.bed}; f=${f%.bim}; f=${f%.fam}; echo "$f"; }
abspath() { local d; d=$(cd "$(dirname "$1")" && pwd) && echo "${d}/$(basename "$1")"; }

file1=$(strip_ext "$1")
file2=$(strip_ext "$2")
file3=$(strip_ext "$3")

mkdir -p "$(dirname "$file3")"
file1=$(abspath "$file1")
file2=$(abspath "$file2")
file3=$(abspath "$file3")

for f in "$file1" "$file2"; do
	for ext in bed bim fam; do
		[ -s "${f}.${ext}" ] || { echo "ERROR: ${f}.${ext} not found" >&2; exit 1; }
	done
done

# All intermediates go in a private temp dir next to the output (not the cwd)
tmp=$(mktemp -d "${file3}.tmp.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

#find common snps
comm -12 <(cut -f 2 ${file1}.bim | sort | uniq) <(cut -f 2 ${file2}.bim | sort | uniq) > ${tmp}/common_snps.txt

#extracting common snps
plink --bfile $file1 --extract ${tmp}/common_snps.txt --allow-no-sex --make-bed --out ${tmp}/temp1_common
plink --bfile $file2 --extract ${tmp}/common_snps.txt --allow-no-sex --make-bed --out ${tmp}/temp2_common

#merging
plink --bfile ${tmp}/temp1_common --bmerge ${tmp}/temp2_common --allow-no-sex --make-bed --out $file3

#missnp file?
if [ -f "${file3}-merge.missnp" ]; then
	echo "Triallelic variants may be present"
	plink --bfile ${tmp}/temp1_common --exclude ${file3}-merge.missnp --make-bed --out ${tmp}/temp1_updated
	plink --bfile ${tmp}/temp2_common --exclude ${file3}-merge.missnp --make-bed --out ${tmp}/temp2_updated

	plink --bfile ${tmp}/temp1_updated --bmerge ${tmp}/temp2_updated.bed ${tmp}/temp2_updated.bim ${tmp}/temp2_updated.fam --allow-no-sex --make-bed --out $file3
fi

echo "done"
