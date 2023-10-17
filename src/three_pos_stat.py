import pandas as pd
import statsmodels.api as sm
from pyfaidx import Fasta
from patsy import dmatrices
import statsmodels.formula.api as smf
import ray

def c_pos(subtype, chromosome, pop = "ALL", suffix = ""):
  """Get the positions for controls for a given subtype"""
  input_dir = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/controls/{}/pos_files/".format(pop)
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt" + suffix
  pos_list = pd.read_csv(f_name, header=None, names = ['pos'], usecols=['pos']).squeeze("columns")
  return pos_list

def s_pos(subtype, chromosome, pop = "ALL"):
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
def get_count_table_control(chromosome, subtype, p, q, r, pop = "ALL"):
  """
  Get the count table for controls at a pair of given relative positions for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
  ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  fasta_obj = Fasta(ref_file)
  control_pos = c_pos(subtype, chromosome, pop)
  rev_control_pos = c_pos(subtype + "_rev", chromosome, pop)
  seq = fasta_obj["{}{}".format("chr", chromosome)]
  seqstr = seq[0:len(seq)].seq
  results = {"AAA":0, "AAC":0, "AAG":0, "AAT":0, 
    "ACA":0, "ACC":0, "ACG":0, "ACT":0,
    "AGA":0, "AGC":0, "AGG":0, "AGT":0,
    "ATA":0, "ATC":0, "ATG":0, "ATT":0,
    "CAA":0, "CAC":0, "CAG":0, "CAT":0, 
    "CCA":0, "CCC":0, "CCG":0, "CCT":0,
    "CGA":0, "CGC":0, "CGG":0, "CGT":0,
    "CTA":0, "CTC":0, "CTG":0, "CTT":0,
    "GAA":0, "GAC":0, "GAG":0, "GAT":0, 
    "GCA":0, "GCC":0, "GCG":0, "GCT":0,
    "GGA":0, "GGC":0, "GGG":0, "GGT":0,
    "GTA":0, "GTC":0, "GTG":0, "GTT":0,
    "TAA":0, "TAC":0, "TAG":0, "TAT":0, 
    "TCA":0, "TCC":0, "TCG":0, "TCT":0,
    "TGA":0, "TGC":0, "TGG":0, "TGT":0,
    "TTA":0, "TTC":0, "TTG":0, "TTT":0
    }
  for index, value in control_pos.items():
    ix1 = value - 1 + p
    ix2 = value - 1 + q
    ix3 = value - 1 + r
    nuc = "{}{}{}".format(seqstr[ix1],seqstr[ix2],seqstr[ix3])
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_control_pos.items():
    ix1 = value - 1 + (p * -1)
    ix2 = value - 1 + (q * -1)
    ix3 = value - 1 + (r * -1)
    nuc = "{}{}{}".format(complement(seqstr[ix1]),complement(seqstr[ix2]),complement(seqstr[ix3]))
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

@ray.remote
def get_count_table_singletons(chromosome, subtype, p, q, r, pop = "ALL"):
  """
  Get the count table for singletons at a given relative position for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
  singleton_pos = s_pos(subtype, chromosome, pop)
  rev_singleton_pos = s_pos(subtype + "_rev", chromosome, pop)
  ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  fasta_obj = Fasta(ref_file)
  seq = fasta_obj["{}{}".format("chr", chromosome)]
  seqstr = seq[0:len(seq)].seq
  results = {"AAA":0, "AAC":0, "AAG":0, "AAT":0, 
    "ACA":0, "ACC":0, "ACG":0, "ACT":0,
    "AGA":0, "AGC":0, "AGG":0, "AGT":0,
    "ATA":0, "ATC":0, "ATG":0, "ATT":0,
    "CAA":0, "CAC":0, "CAG":0, "CAT":0, 
    "CCA":0, "CCC":0, "CCG":0, "CCT":0,
    "CGA":0, "CGC":0, "CGG":0, "CGT":0,
    "CTA":0, "CTC":0, "CTG":0, "CTT":0,
    "GAA":0, "GAC":0, "GAG":0, "GAT":0, 
    "GCA":0, "GCC":0, "GCG":0, "GCT":0,
    "GGA":0, "GGC":0, "GGG":0, "GGT":0,
    "GTA":0, "GTC":0, "GTG":0, "GTT":0,
    "TAA":0, "TAC":0, "TAG":0, "TAT":0, 
    "TCA":0, "TCC":0, "TCG":0, "TCT":0,
    "TGA":0, "TGC":0, "TGG":0, "TGT":0,
    "TTA":0, "TTC":0, "TTG":0, "TTT":0
    }
  for index, value in singleton_pos.items():
    ix1 = value - 1 + p
    ix2 = value - 1 + q
    ix3 = value - 1 + r
    nuc = "{}{}{}".format(seqstr[ix1],seqstr[ix2],seqstr[ix3])
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_singleton_pos.items():
    ix1 = value - 1 + (p * -1)
    ix2 = value - 1 + (q * -1)
    ix3 = value - 1 + (r * -1)
    nuc = "{}{}{}".format(complement(seqstr[ix1]),complement(seqstr[ix2]), complement(seqstr[ix3]))
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

def get_count_all(subtype, p, q, r, status = "singleton", pop = "ALL"):
  """
  Get the singleton and control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
  if status == "singleton":
    futures = [get_count_table_singletons.remote(i, subtype, p, q, r, pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  else:
    futures = [get_count_table_control.remote(i, subtype, p, q, r, pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

def fit_model_all(subtype, p, q, r, pop = "ALL"):
  df_s = get_count_all(subtype, p, q, r, pop = pop)
  df_c = get_count_all(subtype, p, q, r, status = "control", pop = pop)
  df_s = df_s.rename_axis("motif").reset_index()
  df_c = df_c.rename_axis("motif").reset_index()
  df_s.columns = ["motif", "singletons"]
  df_c.columns = ["motif", "controls"]
  df_s['p1'] = df_s['motif'].apply(lambda x: x[0])
  df_s['p2'] = df_s['motif'].apply(lambda x: x[1])
  df_s['p3'] = df_s['motif'].apply(lambda x: x[2])
  df_c['p1'] = df_c['motif'].apply(lambda x: x[0])
  df_c['p2'] = df_c['motif'].apply(lambda x: x[1])
  df_c['p3'] = df_c['motif'].apply(lambda x: x[2])
  df_s = df_s.drop(columns = 'motif')
  df_c = df_c.drop(columns = 'motif')
  df = pd.DataFrame.merge(df_s, df_c, on=['p1', 'p2', 'p3'])
  df = pd.melt(df, id_vars = ['p1', 'p2', 'p3'], var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + p1 + p2 + p1*p2*p3 + p1*p2*status + p2*p3*status + p1*p3*status", df, family = sm.families.Poisson()).fit()
  n_s = sum(df_s.singletons)
  n_c = sum(df_c.controls)
  return {"dev":mod.deviance, "singletons":n_s, "controls":n_c, "rp1":p, "rp2":q, "rp3": r}

ray.init(num_cpus=22)
results = []
pop = "ALL"
subtype = "cpg_GC_CG"
print("Running models for subtype: {} in population: {}".format(subtype, pop))

for p1 in range(-5,4):
    if p1 == 0:
        continue
    if subtype.startswith("cpg") and p1 == 1:
        continue
    print("p1: {}".format(p1))
    for p2 in range((p1+1),5):
        if p2 == 0:
            continue
        if subtype.startswith("cpg") and p2 == 1:
            continue
        print("p2: {}".format(p2))
        for p3 in range((p2+1),6):
            if p3 == 0:
                continue
            if subtype.startswith("cpg") and p3 == 1:
                continue
            print("p3: {}".format(p3))
            results.append(fit_model_all(subtype, p1, p2, p3, pop = pop))

ray.shutdown()
final = pd.DataFrame.from_dict(results)
out_dir = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/three_pos/{}/".format(pop)
file_name = out_dir + subtype + ".csv"
final.to_csv(file_name, index = False)
