from pyfaidx import Fasta
import argparse
import pandas as pd
from Bio.Seq import Seq

def get_motif(seqstr, pos, bp = 1):
    """Get the motif centered at a position from VCF"""
    pyix = pos - 1 # python is index 0, versus 1 for VCF
    return seqstr[(pyix - bp):(pyix + bp + 1)]

def full_cat(ref, alt, motif, bp = 1):
    nuc_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    ref_rc = nuc_dict[ref]
    alt_rc = nuc_dict[alt]
    if ref in ['A','G']:
        final = "{}{}_{}{}".format(ref, ref_rc, alt, alt_rc)
    else:
        final = "{}{}_{}{}".format(ref_rc, ref, alt_rc, alt)
    if motif[bp:bp+2] == "CG":
        final = "cpg_" + final
    return final

def full_motif(motif, ref):
    motif_seq = Seq(motif)
    motif_rc = motif_seq.reverse_complement().__str__()
    if ref in ['A','C']:
        final = "{}({})".format(motif, motif_rc)
    else:
        final = "{}({})".format(motif_rc, motif)
    return final

parser = argparse.ArgumentParser(description="Annotate genomic locations with 21-mer motif")
parser.add_argument("-c", "--chrom", help="Which chromosome?", required=True)
parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
parser.add_argument("-o", "--output", help="Location of output file", required=True)
args = parser.parse_args()

chrom = args.chrom
singleton_file = args.singleton
output_file = args.output

ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
print("Reading reference file...")
fasta_obj = Fasta(ref_file)
seq = fasta_obj["chr{}".format(chrom)]
seqstr = seq[0:len(seq)].seq
print("Reference read!")

output_list = []

with open(singleton_file) as fp:
    line = fp.readline()
    line = fp.readline()
    cnt = 1
    while line:
        data = line.strip().split("\t") # CHROM, POS, S/D, ALLELE, ID
        chrom = data[0]
        pos = int(data[1])
        alt = data[3]
        subject = data[4]
        motif = get_motif(seqstr, pos)
        ref = seqstr[(pos-1)]
        cat = "{}>{}".format(ref, alt)
        motif_full = full_motif(motif, ref)
        entry = {
            'chrom' : chrom,
            'pos' : pos,
            'motif' : motif,
            'cat': cat,
            'alt': alt,
            'subject': subject,
            'ref': ref,
            'motif_full': motif_full
        }
        output_list.append(entry)
        if cnt % 10000 == 0:
            pd.DataFrame(output_list).to_csv(output_file, index = None, header=False, mode='a')
            output_list = []
        line = fp.readline()

if output_list:
    pd.DataFrame(output_list).to_csv(output_file, index = None, header=False, mode='a')
