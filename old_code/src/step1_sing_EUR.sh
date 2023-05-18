#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=700MB
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --job-name=step1_EUR
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step1_eur-%A_%a.err
#SBATCH -o slurm/step1_eur-%A_%a.out

bcftools view --types 'snps' -f 'FILTER="PASS"' -O u /net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/old_code/vcf/1kGP_high_coverage_Illumina.chr${SLURM_ARRAY_TASK_ID}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz |\
bcftools view -S /net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/EUR_IDs.txt -x -Ov |\
vcftools --vcf - --singletons -c > /net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/old_code/output/singletons/EUR/chr${SLURM_ARRAY_TASK_ID}.txt
