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
  seqstr = seq[0:len(seq)].seq
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

ray.init(num_cpus=22)
results = []
pop1 = "EUR"
pop2 = "SAS"
subtype = "cpg_GC_CG"
is_cpg = subtype.startswith("cpg")
print("Running models for subtype: {} in populations: {} and {}".format(subtype, pop1, pop2))

for offset in range(-100, 101):
  if offset == 0:
    continue
  if is_cpg and offset == 1:
    continue
  print(offset)
  df_s_1 = get_count_all(subtype, offset, pop = pop1)
  df_s_2 = get_count_all(subtype, offset, pop = pop2)
  df_c_1 = get_count_all(subtype, offset, status = "control", pop = pop1)
  df_c_2 = get_count_all(subtype, offset, status = "control" , pop = pop2)
  
  df_s_1 = df_s_1.reset_index(level=0)
  df_c_1 = df_c_1.reset_index(level=0)
  df_s_2 = df_s_2.reset_index(level=0)
  df_c_2 = df_c_2.reset_index(level=0)
  
  df_s_1.columns = ['nuc', 'singletons']
  df_s_2.columns = ['nuc', 'singletons']
  df_c_1.columns = ['nuc', 'controls']
  df_c_2.columns = ['nuc', 'controls']
  
  df1 = pd.DataFrame.merge(df_s_1, df_c_1, on='nuc')
  df2 = pd.DataFrame.merge(df_s_2, df_c_2, on='nuc')
  df1 = df1.melt(id_vars = 'nuc', var_name = 'status', value_name = "n")
  df2 = df2.melt(id_vars = 'nuc', var_name = 'status', value_name = "n")
  df1['pop'] = pop1
  df2['pop'] = pop2
  df = pd.concat([df1, df2])
  mod = smf.glm("n ~ status + nuc + pop + status:nuc + pop:nuc", df, family = sm.families.Poisson()).fit()
  
  results.append({"rp": offset, "dev":mod.deviance, "n":df['n'].sum()})

ray.shutdown()
final = pd.DataFrame.from_dict(results)
file_name = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/ref_pop_sp/pop_param/{}_{}_{}.csv".format(subtype, pop1, pop2)
final.to_csv(file_name, index = False)
