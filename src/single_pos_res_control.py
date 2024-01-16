"""
Code for getting the deviance statistics at all positions within +/- 100bp window
"""
import pandas as pd
import statsmodels.api as sm
from pyfaidx import Fasta
from patsy import dmatrices
import statsmodels.formula.api as smf
import ray

def c_pos(subtype, chromosome, pop = "ALL"):
  """Get the positions for controls for a given subtype"""
  input_dir = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/controls/" + pop + "/pos_files/"
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=None, names = ['pos'], usecols=['pos']).squeeze("columns")
  return pos_list

def complement(nucleotide):
  complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
  if not nucleotide in complements.keys():
    return nucleotide
  return complements[nucleotide]

@ray.remote
def get_count_table_control(chromosome, subtype, offset, pop = "ALL"):
  """
  Get the count table for controls at a given relative position for a subtype.
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

def get_count_all(subtype, offset, pop = "ALL"):
  """
  Get the singleton and control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
  futures = [get_count_table_control.remote(i, subtype, offset, pop = pop) for i in range(1,23)]
  results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

def fit_model_all(subtype, offset, pop1, pop2):
  s_tab = get_count_all(subtype, offset, pop = pop1)
  s_tab = s_tab.reset_index(level=0)
  s_tab.columns = ['nuc', 'singletons']
  c_tab = get_count_all(subtype, offset, pop = pop2)
  c_tab = c_tab.reset_index(level=0)
  c_tab.columns = ['nuc', 'controls']
  df = pd.DataFrame.merge(s_tab, c_tab, on='nuc')
  df = pd.melt(df, id_vars = 'nuc', var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + nuc", df, family = sm.families.Poisson()).fit()
  df['fitted'] = mod.fittedvalues
  df['res'] = mod.resid_deviance
  return(df)

ray.init(num_cpus=22)
subtype = "AT_TA"
pop1 = "AFR"
pop2 = "EUR"
out_dir = "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/single_pos/resid/control_control/"
print("Running models for subtype: {} and populations: {} and {}".format(subtype, pop1, pop2))
for offset in range(1, 101):
  print(offset)
  df = fit_model_all(subtype, offset * -1, pop1 = pop1, pop2 = pop2)
  file_name = out_dir + subtype + "_" + pop1 + "_" + pop2 + "_rp_" + str(-1*offset) + ".csv"
  df.to_csv(file_name, index = False)
  if subtype.startswith("cpg") and offset == 1:
    continue
  df = fit_model_all(subtype, offset, pop1 = pop1, pop2 = pop2)
  file_name = out_dir + subtype + "_" + pop1 + "_" + pop2 + "_rp_" + str(offset) + ".csv"
  df.to_csv(file_name, index = False)

ray.shutdown()


