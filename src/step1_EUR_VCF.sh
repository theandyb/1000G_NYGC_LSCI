#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=700MB
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --array=1-22
#SBATCH --job-name=singEUR
#SBATCH --requeue
#SBATCH -e slurm/singleton_eurVCF-%J.err
#SBATCH -o slurm/singleton_eurVCF-%J.out

VCF_DIR="/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/data"
OUT_DIR="/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/data"

bcftools view -i "%FILTER=='PASS'" -v snps ${VCF_DIR}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${SLURM_ARRAY_TASK_ID}.recalibrated_variants.annotated.vcf.gz |\
bcftools view -i "(AC_EUR_unrel == 1 | (AC_Hom_EUR_unrel == 2 & AC_Het_EUR_unrel ==0)) & AC_AFR_unrel == 0 & AC_AMR_unrel == 0 & AC_EAS_unrel == 0 & AC_SAS_unrel == 0" >\
$OUT_DIR/EUR_${SLURM_ARRAY_TASK_ID}_singletons.vcf
