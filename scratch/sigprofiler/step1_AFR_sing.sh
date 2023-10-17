#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=700MB
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --job-name=singAFR
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/singleton_afr-%J.err
#SBATCH -o slurm/singleton_afr-%J.out

VCF_DIR="/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/data/"
OUT_DIR="/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/sigprofiler/output/vcf/"

bcftools view -v snps ${VCF_DIR}/1kGP_high_coverage_Illumina.chr${SLURM_ARRAY_TASK_ID}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz |\
bcftools view -S $VCF_DIR/AFR_IDs.txt |\
bcftools annotate -x INFO,^FORMAT/GT |\
vcftools --singletons --vcf - --out chr${SLURM_ARRAY_TASK_ID}_afr_sing

tail -n +2 chr${SLURM_ARRAY_TASK_ID}_afr_sing.singletons | awk '{print($1"\t"$2)}' > sites_afr_${SLURM_ARRAY_TASK_ID}.txt

bcftools view -v snps ${VCF_DIR}/1kGP_high_coverage_Illumina.chr${SLURM_ARRAY_TASK_ID}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz |\
bcftools view -S $VCF_DIR/AFR_IDs.txt -T sites_afr_${SLURM_ARRAY_TASK_ID}.txt  |\
bcftools annotate -Ob -x INFO,^FORMAT/GT > $OUT_DIR/afr_chr_${SLURM_ARRAY_TASK_ID}.bcf
