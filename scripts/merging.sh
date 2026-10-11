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

file1=$1
file2=$2
file3=$3
#find common snps
#module load python/3.8.6 R
#Rscript commonsnps.R ${file1}.bim ${file2}.bim

comm -12 <(cut -f 2 ${file1}.bim | sort | uniq) <(cut -f 2 ${file2}.bim | sort | uniq) > common_snps.txt

#extracting common snps
plink --bfile $file1 --extract common_snps.txt --allow-no-sex --make-bed --out temp1_common
plink --bfile $file2 --extract common_snps.txt --allow-no-sex  --make-bed --out temp2_common

#merging
plink --bfile temp1_common --bmerge temp2_common --allow-no-sex --make-bed --out $file3

#missnp file?
if [ -f "${file3}-merge.missnp" ]; then
	echo "Triallelic variants may be present"
	plink --bfile temp1_common --exclude ${file3}-merge.missnp --make-bed --out temp1_updated
	plink --bfile temp2_common --exclude ${file3}-merge.missnp --make-bed --out temp2_updated

	plink --bfile temp1_updated --bmerge temp2_updated.bed temp2_updated.bim temp2_updated.fam --allow-no-sex --make-bed --out $file3
fi

echo "done"
rm temp*
