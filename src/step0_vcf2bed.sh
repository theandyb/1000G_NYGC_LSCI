#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=2GB
#SBATCH --ntasks=1
#SBATCH --time 02:00:00
#SBATCH --job-name=step0
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/vcf_loc-%J.err
#SBATCH -o slurm/vcf_loc-%J.out

VCF_DIR="/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/data"

bcftools view -G --types 'snps' -f 'FILTER="PASS"' -Ov ${VCF_DIR}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${SLURM_ARRAY_TASK_ID}.recalibrated_variants.annotated.vcf.gz |\
bcftools query -f '%CHROM\t%POS\t%POS\n' | awk '{print($1"\t"$2 - 1"\t"$3)}' >\
/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/data/vcfbed/chr${SLURM_ARRAY_TASK_ID}.bed

