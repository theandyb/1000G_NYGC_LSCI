"""
Code for getting the deviance statistics at all positions within +/- 1000bp window
"""
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
def get_count_table_control(chromosome, subtype, offset, pop = "ALL", suffix = ""):
  """
  Get the count table for controls at a given relative position for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
  ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  fasta_obj = Fasta(ref_file)
  control_pos = c_pos(subtype, chromosome, pop, suffix = suffix)
  rev_control_pos = c_pos(subtype + "_rev", chromosome, pop, suffix = suffix)
  seq = fasta_obj["{}{}".format("chr", chromosome)]
  # seqstr = seq[0:len(seq)].seq
  results = {"A":0, "C":0, "G":0, "T":0}
  for index, value in control_pos.items():
    ix = value - 1 + offset
    nuc = seqstr[ix]
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_control_pos.items():
    ix = value - 1 + (offset * -1)
    nuc = complement(seqstr[ix])
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

@ray.remote
def get_count_table_singletons(chromosome, subtype, offset, pop = "ALL"):
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

def get_count_all(subtype, offset, status = "singleton", pop = "ALL", suffix = ""):
  """
  Get the singleton and control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
  if status == "singleton":
    futures = [get_count_table_singletons.remote(i, subtype, offset, pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  else:
    futures = [get_count_table_control.remote(i, subtype, offset, pop, suffix = suffix) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

def fit_model_all(subtype, offset, pop = "ALL", suffix = ""):
  s_tab = get_count_all(subtype, offset, status = "singleton", pop = pop, suffix = suffix)
  s_tab = s_tab.reset_index(level=0)
  s_tab.columns = ['nuc', 'singletons']
  c_tab = get_count_all(subtype, offset, status = "control", pop = pop, suffix = suffix)
  c_tab = c_tab.reset_index(level=0)
  c_tab.columns = ['nuc', 'controls']
  df2 = pd.DataFrame.merge(s_tab, c_tab, on='nuc')
  df = pd.melt(df2, id_vars = 'nuc', var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + nuc", df, family = sm.families.Poisson()).fit()
  n_s = sum(s_tab.singletons)
  n_c = sum(c_tab.controls)
  df2['p'] = df2['singletons'] / df2['singletons'].sum()
  df2['q'] = df2['controls'] / df2['controls'].sum()
  df2['tv'] = abs(df2['p'] - df2['q'])
  tot_var = df2['tv'].sum()
  max_d = round(df2['tv'].max(),5)
  mean_var = df2['tv'].mean()
  return {"dev":mod.deviance, "singletons":n_s, "controls":n_c, "offset":offset, "tot_var":tot_var,"mean_var": mean_var ,"max_var": max_d}

ray.init(num_cpus=22)
results = []
pop = "ALL"
subtype = "GC_TA"
suffix = ".max"
print("Running models for subtype: {} in population: {}".format(subtype, pop))
for offset in range(1, 1001):
  print(offset)
  results.append(fit_model_all(subtype, offset * -1, pop = pop, suffix = suffix))
  if subtype.startswith("cpg") and offset == 1:
    continue
  results.append(fit_model_all(subtype, offset, pop = pop, suffix = suffix))

ray.shutdown()
  
final = pd.DataFrame.from_dict(results)
out_dir = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/single_pos/{}/".format(pop)
file_name = out_dir + subtype + ".csv" + suffix
final.to_csv(file_name, index = False)

