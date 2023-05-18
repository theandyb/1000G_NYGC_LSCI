from pyfaidx import Fasta
import argparse
import pandas as pd
#from Bio.Seq import Seq

parser = argparse.ArgumentParser(description="Identify CpG and Non-CpG singletons")
parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
parser.add_argument("-o", "--output", help="Location of output file", required=True)
args = parser.parse_args()

singleton_file = args.singleton
output_file = args.output

ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
fasta_obj = Fasta(ref_file)

current_chrom = 1
seq = fasta_obj["chr{}".format(current_chrom)]
seqstr = seq[0:len(seq)].seq

output_list = []

with open(singleton_file) as fp:
  cnt = 1
  line = fp.readline()
  while(line):
    data = line.strip().split("\t") # CHROM, POS, REF
    chrom = int(data[0][3:])
    pos = int(data[1])
    ref_al = data[2]
    if(current_chrom != chrom):
      current_chrom = chrom
      seq = fasta_obj["chr{}".format(current_chrom)]
      seqstr = seq[0:len(seq)].seq
    ref = seqstr[(pos-1):(pos+1)] # G upstream of a C?
    ref2 = seqstr[(pos-2):(pos)] # C downstream of a G?
    if(ref == "CG" or ref2 == "CG"):
      cpg_stat = 1
    else:
      cpg_stat = 0
    entry = {'chrom': chrom, 'pos': pos, 'ref': ref_al, 'cpg': cpg_stat}
    output_list.append(entry)
    cnt += 1
    if cnt % 10000 == 0:
      pd.DataFrame(output_list).to_csv(output_file, index = None, header=False, mode='a')
      output_list = []
    line = fp.readline()

if output_list:
    pd.DataFrame(output_list).to_csv(output_file, index = None, header=False, mode='a')
