import pandas as pd
from pyfaidx import Fasta
import ray

def s_pos(subtype, chromosome, pop):
    """Get the positions for singletons for a given subtype"""
    input_dir = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/{}/pos_files/".format(pop)
    f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
    pos_list = pd.read_csv(f_name, header=None, names = ['pos'], usecols=['pos']).squeeze("columns")
    return pos_list

def complement(nucleotide):
    complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    if not nucleotide in complements.keys():
        return nucleotide
    return complements[nucleotide]

@ray.remote
def get_count_table(chromosome, subtype, offset, pop):
    singleton_pos = s_pos(subtype, chromosome, pop)
    rev_singleton_pos = s_pos(subtype + "_rev", chromosome, pop)
    ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    fasta_obj = Fasta(ref_file)
    seq = fasta_obj["{}{}".format("chr", chromosome)]
    seqstr = seq[0:len(seq)].seq
    results = {"A":0, "C":0, "G":0, "T":0}
    for index, value in singleton_pos.items():
        ix = value - 1 + offset
        nuc = seqstr[ix]
        if nuc in results.keys():
            results[nuc] += 1
    for index, value in rev_singleton_pos.items():
        ix = value - 1 + (offset * -1)
        nuc = complement(seqstr[ix])
        if nuc in results.keys():
            results[nuc] += 1
    return(results)

def get_count_all(subtype, offset, pop):
    futures = [get_count_table.remote(i, subtype, offset, pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
    return(results)

def get_stat(pop1, pop2, subtype, offset):
    df1 = get_count_all(subtype, offset, pop1)
    df2 = get_count_all(subtype, offset, pop2)
    p_1 = df1 / sum(df1)
    e_2 = p_1 * sum(df2)
    chi_2 = ((df2 - e_2)**2) / e_2
    stat_val = sum(chi_2)
    return(stat_val)

ray.init(num_cpus=22)

pop1 = "SAS" # reference population
pop2 = "EUR"

results = []
subtype = "cpg_GC_CG"
print("Computing statistics for {} using {} as reference population, subtype {}".format(pop2, pop1, subtype))
for offset in range(1,101):
    print(offset)
    results.append({"rp":(-1*offset), "stat": get_stat(pop1, pop2, subtype, (-1 * offset))})
    if subtype.startswith("cpg") and offset == 1:
        continue
    results.append({"rp":offset, "stat": get_stat(pop1, pop2, subtype, offset)})

ray.shutdown()

final = pd.DataFrame.from_dict(results)
out_file = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/ref_pop_sp/{}_{}_{}.csv".format(pop1, pop2, subtype)
final.to_csv(out_file, index = False)
