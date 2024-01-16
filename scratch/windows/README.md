# Window Analysis

## Introduction

The goal of this analysis is to assess the consistency of our single position results in windows across the genome. We'll first consider 50kb windows, as this was the size utilized by Seplyarskiy, et al (2023) in defining local variation in mutation rate.

## First pass

We'll first look at the singletons (and their matched controls) from the EUR superpopulation for the A to G mutation subtype. To make things simpler (to begin), we'll first look at chromosome 22.

### Get singletons in 50kb windows

First lets assign each singleton of type A to G in chr22 to a bin and see what the distribution of counts per bin looks like:

```
awk 'function ceil(x, y){y=int(x); return(x>y?y+1:y)} {if($1=="chr22")count[ceil($2/50000)]++}END{for(key in count)print(key"\t"count[key])}' ../../output/singletons/EUR/AT_GC.txt > bin_counts.txt
```

Over 96% of bins have at least 10 singleton observations, so we will proceed with performing our models using 50kb as our bin size.

```
awk 'function ceil(x, y){y=int(x); return(x>y?y+1:y)} {if($1=="chr22")print($1"\t"$2"\t"$3"\t"ceil($2/50000))}' ../../output/singletons/EUR/AT_GC.txt > data/singletons/AT_GC_EUR_22.txt
```

#### Position files

```
awk '{if($3=="A")print($2"\t"$4)}' AT_GC_EUR_22.txt > pos_files/AT_GC_EUR_22.txt
awk '{if($3=="T")print($2"\t"$4)}' AT_GC_EUR_22.txt > pos_files/AT_GC_EUR_22_rev.txt
```

### Get controls in 50kb windows

We will define windows based on where matched singletons are located

```
awk -F, 'function ceil(x, y){y=int(x); return(x>y?y+1:y)} {if($1 == "chr22") print($1"\t"$7"\t"$3"\t"ceil($2/50000)"\t"$6"\t"$2)}' ../../output/controls/EUR/AT_GC.csv > data/controls/AT_GC_EUR_22.txt
```

#### Position files

```
awk '{if($3=="A")print($2"\t"$4"\t"$5"\t"$6)}' AT_GC_EUR_22.txt > pos_files/AT_GC_EUR_22.txt
awk '{if($3=="T")print($2"\t"$4"\t"$5"\t"$6)}' AT_GC_EUR_22.txt > pos_files/AT_GC_EUR_22_rev.txt
```




