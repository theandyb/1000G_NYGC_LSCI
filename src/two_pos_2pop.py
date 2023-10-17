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

def fit_model_all(subtype, p, q, pop1, pop2):
  df_s1 = get_count_all(subtype, p, q, pop = pop1)
  df_s2 = get_count_all(subtype, p, q, pop = pop2)
  df_c1 = get_count_all(subtype, p, q, status = "control", pop = pop1)
  df_c2 = get_count_all(subtype, p, q, status = "control", pop = pop2)
  df_s1 = df_s1.rename_axis("motif").reset_index()
  df_s2 = df_s2.rename_axis("motif").reset_index()
  df_c1 = df_c1.rename_axis("motif").reset_index()
  df_c2 = df_c2.rename_axis("motif").reset_index()
  df_s1.columns = ["motif", "singletons"]
  df_s2.columns = ["motif", "singletons"]
  df_c1.columns = ["motif", "controls"]
  df_c2.columns = ["motif", "controls"]
  df_s1['p1'] = df_s1['motif'].apply(lambda x: x[0])
  df_s2['p1'] = df_s2['motif'].apply(lambda x: x[0])
  df_s1['p2'] = df_s1['motif'].apply(lambda x: x[1])
  df_s2['p2'] = df_s2['motif'].apply(lambda x: x[1])
  df_c1['p1'] = df_c1['motif'].apply(lambda x: x[0])
  df_c2['p1'] = df_c2['motif'].apply(lambda x: x[0])
  df_c1['p2'] = df_c1['motif'].apply(lambda x: x[1])
  df_c2['p2'] = df_c2['motif'].apply(lambda x: x[1])
  df_s1 = df_s1.drop(columns = 'motif')
  df_s2 = df_s2.drop(columns = 'motif')
  df_c1 = df_c1.drop(columns = 'motif')
  df_c2 = df_c2.drop(columns = 'motif')
  df_s1['pop'] = pop1
  df_s2['pop'] = pop2
  df_c1['pop'] = pop1
  df_c2['pop'] = pop2
  df1 = pd.DataFrame.merge(df_s1, df_c1, on=['p1', 'p2', 'pop'])
  df2 = pd.DataFrame.merge(df_s2, df_c2, on=['p1', 'p2', 'pop'])
  df = pd.concat([df1, df2])
  df = pd.melt(df, id_vars = ['p1', 'p2', 'pop'], var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ (p1 * p2 * status) + pop + pop:p1 + pop:p2 + pop:p1:status + pop:p2:status", df, family = sm.families.Poisson()).fit()
  n_s = sum(df_s1.singletons) + sum(df_s2.singletons)
  n_c = sum(df_c1.controls) + sum(df_c2.controls)
  return {"dev":mod.deviance, "singletons":n_s, "controls":n_c, "rp1":p, "rp2":q}

ray.init(num_cpus=22)
results = []
pop1 = "EAS"
pop2 = "EUR"
subtype = "cpg_GC_TA"
print("Running models for subtype: {} in populations: {} and {}".format(subtype, pop1, pop2))

for p1 in range(-10, 10):
  print("p1: {}".format(p1))
  if p1 == 0:
    continue
  if subtype.startswith("cpg") and p1 == 1:
    continue
  for p2 in range((p1+1),11):
    if p2 == 0: 
      continue
    if subtype.startswith("cpg") and p2 == 1:
      continue
    print("p2: {}".format(p2))
    results.append(fit_model_all(subtype, p1, p2, pop1, pop2))

ray.shutdown()
final = pd.DataFrame.from_dict(results)
out_dir = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/two_pos_2pop/"
file_name = out_dir + pop1 + "_" + pop2 + "_" + subtype + ".csv"
final.to_csv(file_name, index = False)
