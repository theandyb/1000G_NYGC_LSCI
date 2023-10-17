import pandas as pd
from pyfaidx import Fasta
from Bio.Seq import Seq
import itertools

ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
print("Reading reference file...")
fasta_obj = Fasta(ref_file)
# chrom = "chr1"
# current_chrom = chrom
# seq = fasta_obj[chrom]
# seqstr = seq[0:len(seq)].seq

nucs = ["A", "C", "G", "T"]

motifs = [''.join(i) for i in list(itertools.product(*[nucs, nucs, nucs]))]
results = {key:0 for key in motifs}

for chrom in range(1, 23):
    print(chrom)
    chrom_name = "chr{}".format(chrom)
    seq = fasta_obj[chrom_name]
    seqstr = seq[0:len(seq)].seq
    for i in range(len(seq)-2):
        test_motif = seqstr[i:(i+3)]
        if test_motif in motifs:
            results[test_motif] += 1

df = pd.DataFrame.from_dict(results, 'index', columns = ["count"])
df.index.name = 'motif'
df.reset_index(inplace = True)
df.to_csv("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/gw_3mer.csv")
