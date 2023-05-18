# get the data
for i in `seq 2 22`; do
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.recalibrated_variants.annotated.vcf.gz
end

# Attempt to get singleton (and doubleton) locations from VCF for SAS
bcftools view -i "%FILTER=='PASS'" -v snps 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.recalibrated_variants.annotated.vcf.gz |\
bcftools view -i "(AC_SAS_unrel == 1 | (AC_Hom_SAS_unrel == 2 & AC_Het_SAS_unrel ==0)) & AC_AFR_unrel == 0 & AC_EUR_unrel == 0 & AC_EAS_unrel == 0 & AC_AMR_unrel == 0" |\
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" > /net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/SAS/chr1.txt

# Subtype specific files
for i in `seq 1 22`; do
echo $i
awk '{if(($3 == "A" && $4 == "C")|| ($3 == "T" && $4 == "G"))print($1"\t"$2"\t"$3)}' chr$i.txt >> AT_CG.txt
awk '{if(($3 == "A" && $4 == "G")|| ($3 == "T" && $4 == "C"))print($1"\t"$2"\t"$3)}' chr$i.txt >> AT_GC.txt
awk '{if(($3 == "A" && $4 == "T")|| ($3 == "T" && $4 == "A"))print($1"\t"$2"\t"$3)}' chr$i.txt >> AT_TA.txt
awk '{if(($3 == "C" && $4 == "A")|| ($3 == "G" && $4 == "T"))print($1"\t"$2"\t"$3)}' chr$i.txt >> GC_TA.txt
awk '{if(($3 == "C" && $4 == "G")|| ($3 == "G" && $4 == "C"))print($1"\t"$2"\t"$3)}' chr$i.txt >> GC_CG.txt
awk '{if(($3 == "C" && $4 == "T")|| ($3 == "G" && $4 == "A"))print($1"\t"$2"\t"$3)}' chr$i.txt >> GC_AT.txt
done

# Annotate GC subtypes with CpG status
python src/cpg_status.py -s output/singletons/SAS/GC_AT.txt -o output/singletons/SAS/GC_AT_ann.txt
python src/cpg_status.py -s output/singletons/SAS/GC_TA.txt -o output/singletons/SAS/GC_TA_ann.txt
python src/cpg_status.py -s output/singletons/SAS/GC_CG.txt -o output/singletons/SAS/GC_CG_ann.txt

awk -F, '{if($4==1)print("chr"$1"\t"$2"\t"$3)}' GC_AT_ann.txt > cpg_GC_AT.txt
awk -F, '{if($4==0)print("chr"$1"\t"$2"\t"$3)}' GC_AT_ann.txt > GC_AT.txt
awk -F, '{if($4==1)print("chr"$1"\t"$2"\t"$3)}' GC_TA_ann.txt > cpg_GC_TA.txt
awk -F, '{if($4==0)print("chr"$1"\t"$2"\t"$3)}' GC_TA_ann.txt > GC_TA.txt
awk -F, '{if($4==1)print("chr"$1"\t"$2"\t"$3)}' GC_CG_ann.txt > cpg_GC_CG.txt
awk -F, '{if($4==0)print("chr"$1"\t"$2"\t"$3)}' GC_CG_ann.txt > GC_CG.txt

# Sample controls
python sample_control.py -s ../output/singletons/SAS/AT_CG.txt -f ../data/mask_ref.fa -o ../output/controls/SAS/AT_CG.csv -t "AT_CG" -n 5
python sample_control.py -s ../output/singletons/SAS/AT_GC.txt -f ../data/mask_ref.fa -o ../output/controls/SAS/AT_GC.csv -t "AT_GC" -n 5
python sample_control.py -s ../output/singletons/SAS/AT_TA.txt -f ../data/mask_ref.fa -o ../output/controls/SAS/AT_TA.csv -t "AT_TA" -n 5

python sample_control.py -s ../output/singletons/SAS/GC_AT.txt -f ../data/mask_ref.fa -o ../output/controls/SAS/GC_AT.csv -t "GC_AT" -n 5
python sample_control.py -s ../output/singletons/SAS/GC_TA.txt -f ../data/mask_ref.fa -o ../output/controls/SAS/GC_TA.csv -t "GC_TA" -n 5
python sample_control.py -s ../output/singletons/SAS/GC_CG.txt -f ../data/mask_ref.fa -o ../output/controls/SAS/GC_CG.csv -t "GC_CG" -n 5

python sample_control.py -s ../output/singletons/SAS/cpg_GC_AT.txt -f ../data/mask_ref.fa -o ../output/controls/SAS/cpg_GC_AT.csv -t "cpg_GC_AT" -n 5
python sample_control.py -s ../output/singletons/SAS/cpg_GC_TA.txt -f ../data/mask_ref.fa -o ../output/controls/SAS/cpg_GC_TA.csv -t "cpg_GC_TA" -n 5
python sample_control.py -s ../output/singletons/AMR/cpg_GC_CG.txt -f ../data/mask_ref.fa -o ../output/controls/AMR/cpg_GC_CG.csv -t "cpg_GC_CG" -n 5

# ALL file generation

#################################################
## Paste singleton files into one
sub=""
cat AFR/${sub}.txt >> ${sub}.txt
cat AMR/${sub}.txt >> ${sub}.txt
cat EAS/${sub}.txt >> ${sub}.txt
cat SAS/${sub}.txt >> ${sub}.txt
cat EUR/${sub}.txt >> ${sub}.txt

## Ensure proper sorting
sort -k1,1V -k2,2n ${sub}.txt > ALL/${sub}.txt
#################################################

#################################################
## Paste control files into one
sub="cpg_GC_CG"
cat AFR/${sub}.csv >> ${sub}.csv
cat AMR/${sub}.csv >> ${sub}.csv
cat EAS/${sub}.csv >> ${sub}.csv
cat EUR/${sub}.csv >> ${sub}.csv
cat SAS/${sub}.csv >> ${sub}.csv

## Ensure proper sorting
sort -k1,1V -k2,2n ${sub}.csv > ALL/${sub}.csv
#################################################

# pos_files
## Singletons
sub="AT_TA"
for i in `seq 1 22`; do
echo $i
awk '{if($1=="chr'$i'" && $3 == "A")print($2)}' $sub.txt > pos_files/${sub}_${i}.txt
awk '{if($1=="chr'$i'" && $3 == "T")print($2)}' $sub.txt > pos_files/${sub}_rev_${i}.txt
done

# GC version
sub="GC_AT"
for i in `seq 1 22`; do
echo $i
awk '{if($1=="chr'$i'" && $3 == "C")print($2)}' $sub.txt > pos_files/${sub}_${i}.txt
awk '{if($1=="chr'$i'" && $3 == "G")print($2)}' $sub.txt > pos_files/${sub}_rev_${i}.txt
done


## controls
sub="AT_TA"
for i in `seq 1 22`; do
echo $i
awk -F, '{if($1=="chr'$i'" && $4 == "A")print($7)}' $sub.csv > pos_files/${sub}_${i}.txt
awk -F, '{if($1=="chr'$i'" && $4 == "T")print($7)}' $sub.csv > pos_files/${sub}_rev_${i}.txt
done
# GC version
sub="cpg_GC_CG"
for i in `seq 1 22`; do
echo $i
awk -F, '{if($1=="chr'$i'" && $4 == "C")print($7)}' $sub.csv > pos_files/${sub}_${i}.txt
awk -F, '{if($1=="chr'$i'" && $4 == "G")print($7)}' $sub.csv > pos_files/${sub}_rev_${i}.txt
done


