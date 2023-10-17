"""
Code for getting the deviance statistics for two position models
"""
import pandas as pd
import statsmodels.api as sm
from pyfaidx import Fasta
from patsy import dmatrices
import statsmodels.formula.api as smf
import ray

def c_pos(subtype, chromosome, pop="ALL"):
  """Get the positions for controls for a given subtype"""
  input_dir = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/controls/{}/pos_files/".format(pop)
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=0, names = ['pos'], usecols=['pos']).squeeze("columns")
  return pos_list

def s_pos(subtype, chromosome, pop="ALL"):
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
def get_count_table_control(chromosome, subtype, p, q, pop = "ALL"):
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
  results = {"AA":0, "AC":0, "AG":0, "AT":0, 
    "CA":0, "CC":0, "CG":0, "CT":0,
    "GA":0, "GC":0, "GG":0, "GT":0,
    "TA":0, "TC":0, "TG":0, "TT":0}
  for index, value in control_pos.items():
    ix1 = value - 1 + p
    ix2 = value - 1 + q
    nuc = "{}{}".format(seqstr[ix1],seqstr[ix2])
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_control_pos.items():
    ix1 = value - 1 + (p * -1)
    ix2 = value - 1 + (q * -1)
    nuc = "{}{}".format(complement(seqstr[ix1]),complement(seqstr[ix2]))
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

@ray.remote
def get_count_table_singletons(chromosome, subtype, p, q, pop = "ALL"):
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
  results = {"AA":0, "AC":0, "AG":0, "AT":0, 
    "CA":0, "CC":0, "CG":0, "CT":0,
    "GA":0, "GC":0, "GG":0, "GT":0,
    "TA":0, "TC":0, "TG":0, "TT":0}
  for index, value in singleton_pos.items():
    ix1 = value - 1 + p
    ix2 = value - 1 + q
    nuc = "{}{}".format(seqstr[ix1],seqstr[ix2])
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_singleton_pos.items():
    ix1 = value - 1 + (p * -1)
    ix2 = value - 1 + (q * -1)
    nuc = "{}{}".format(complement(seqstr[ix1]),complement(seqstr[ix2]))
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

def get_count_all(subtype, p, q, status = "singleton", pop = "ALL"):
  """
  Get the singleton and control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
  if status == "singleton":
    futures = [get_count_table_singletons.remote(i, subtype, p, q, pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  else:
    futures = [get_count_table_control.remote(i, subtype, p, q, pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

def fit_model_all(subtype, p, q, pop = "ALL"):
  df_s = get_count_all(subtype, p, q, pop = pop)
  df_c = get_count_all(subtype, p, q, status = "control", pop = pop)
  df_s = df_s.rename_axis("motif").reset_index()
  df_c = df_c.rename_axis("motif").reset_index()
  df_s.columns = ["motif", "singletons"]
  df_c.columns = ["motif", "controls"]
  df_s['p1'] = df_s['motif'].apply(lambda x: x[0])
  df_s['p2'] = df_s['motif'].apply(lambda x: x[1])
  df_c['p1'] = df_c['motif'].apply(lambda x: x[0])
  df_c['p2'] = df_c['motif'].apply(lambda x: x[1])
  df_s = df_s.drop(columns = 'motif')
  df_c = df_c.drop(columns = 'motif')
  df = pd.DataFrame.merge(df_s, df_c, on=['p1', 'p2'])
  df = pd.melt(df, id_vars = ['p1', 'p2'], var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + p1 + p2 + p1*p2 + p1*status + p2*status", df, family = sm.families.Poisson()).fit()
  df['fitted'] = mod.fittedvalues
  df['res'] = mod.resid_deviance
  return df

ray.init(num_cpus=22)
pop = "EUR"
subtype = "GC_AT"
out_dir = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/two_pos/resid/{}/".format(pop)

print("Running models for subtype: {} in population: {}".format(subtype, pop))

for p1 in range(-20, 20):
  print("p1: {}".format(p1))
  if p1 == 0:
    continue
  if subtype.startswith("cpg") and p1 == 1:
    continue
  for p2 in range((p1+1),21):
    if p2 == 0: 
      continue
    if subtype.startswith("cpg") and p2 == 1:
      continue
    print("p2: {}".format(p2))
    df = fit_model_all(subtype, p1, p2, pop = pop)
    file_name = out_dir + subtype + "_p" + str(p1) + "_q" + str(p2) + ".csv"
    df.to_csv(file_name, index = False)

ray.shutdown()

