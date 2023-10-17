cd /net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/data
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz

bcftools view -S unrel_ids.txt -v snps 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz |\
vcftools --vcf - --singletons --out chr22_sing

python ../src/annotate_vcftools_singleton.py -s chr22_sing.singletons -c 22 -o chr22_sing.ann.csv
