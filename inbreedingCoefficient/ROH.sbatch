#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=96:00:00
#SBATCH --job-name=R
#SBATCH --output=R_%j.out 
#SBATCH --error=R_%j.err
#SBATCH --mail-user=lilepisorus@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

cd $SLURM_SUBMIT_DIR
module load vcftools
module load plink

vcftools --vcf allSNPs.vcf --thin 2000 --recode --out maize.thin
#convert to the plink format
vcftools --vcf maize.thin.recode.vcf --plink --out maize.thin

#Run the ROHs analysis
plink --file maize.thin --out ROH300kb.20snpWindow.1het --homozyg-window-snp 20 --homozyg-snp 10 --homozyg-kb 300 --homozyg-gap 1000 --homozyg-window-missing 2 --homozyg-window-threshold .05 --homozyg-window-het 1 --homozyg-density 50

