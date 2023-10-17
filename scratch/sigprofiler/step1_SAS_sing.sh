#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=700MB
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --job-name=singSAS
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/singleton_sas-%J.err
#SBATCH -o slurm/singleton_sas-%J.out

VCF_DIR="/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/data/"
OUT_DIR="/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/SAS/"

bcftools view -i "%FILTER=='PASS'" -v snps ${VCF_DIR}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${SLURM_ARRAY_TASK_ID}.recalibrated_variants.annotated.vcf.gz |\
bcftools view -i "(AC_SAS_unrel == 1 | (AC_Hom_SAS_unrel == 2 & AC_Het_SAS_unrel ==0)) & AC_AFR_unrel == 0 & AC_AMR_unrel == 0 & AC_EAS_unrel == 0 & AC_EUR_unrel == 0" |\
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" > $OUT_DIR/chr${SLURM_ARRAY_TASK_ID}.txt
