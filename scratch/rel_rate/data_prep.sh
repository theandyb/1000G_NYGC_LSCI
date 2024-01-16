# Get 3-mer counts from 1kGP
## singletons
### central
pop="SAS"
awk -F, '{count[substr($4,10,3)]++}END{for(key in count)print(key"\t"count[key])}' /net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/${pop}/motifs/AT_GC.csv >\
/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/1kgp_rates/${pop}/A_G_3mer_count.txt
### augmented
awk -F, '{count[substr($4,9,3)]++}END{for(key in count)print(key"\t"count[key])}' /net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/${pop}/motifs/AT_GC.csv >\
/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/1kgp_rates/${pop}/A_G_3mer_count_aug.txt

## controls
pop="ALL"
awk -F, '{count[substr($4,10,3)]++}END{for(key in count)print(key"\t"count[key])}' /net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/controls/${pop}/motifs/AT_GC.csv >\
/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/1kgp_rates/${pop}/A_G_control_3mer_count.txt
### augmented
awk -F, '{count[substr($4,9,3)]++}END{for(key in count)print(key"\t"count[key])}' /net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/controls/${pop}/motifs/AT_GC.csv >\
/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/1kgp_rates/${pop}/A_G_control_3mer_count_aug.txt

# HGDP singletons

## batch script get_singletons.sh runs vcftools on (slightly) filtered vcfs from hgdp

## annotation script

## hgdp_singleton_files.R to generate population-specific singleton list

# sample hgdp controls
python sample_hgdp_control.py -s hgdp_singletons/Africa_A_G.csv -f ../../data/mask_ref.fa -o hgdp_controls/Africa_A_G.csv -t "AT_GC" -n 2
# annotate
python annotate_control_hgdp.py -c hgdp_controls/Middle_East_A_G.csv -o hgdp_controls/annotated/Middle_East_A_G.csv
