import logging as log
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

sys.path.append("../PrEWInputAnalysis")
import FuncHelp.Wrappers as FHW
import IO.FilenameHelp as IOFH
import IO.Reader as IOR
import IO.SysHelpers as IOSH
import Plotting.DefaultFormat as PDF
import Plotting.Naming as PN
import Shape.ShapeFunctions as SSF
import Shape.ShapeTesting as SST


def dif_param_LR(cos_th, xs0, Ae, Af, ef, k0, dk):
  """ Difermion parametrisation in the difermion rest frame for the left-handed
      electron right-handed positron initial state.
  """
  return 3./8. * xs0 * (1.0 + Ae)/2.0 * ( (1. + (k0 + dk)/2.0) + (ef + 2.0 * Af) * cos_th + (1.0 - 3.0 * (k0 + dk)/2.0) * cos_th*cos_th )

def dif_param_RL(cos_th, xs0, Ae, Af, ef, k0, dk):
  """ Difermion parametrisation in the difermion rest frame for the right-handed
      electron left-handed positron initial state.
  """
  return 3./8. * xs0 * (1.0 - Ae)/2.0 * ( (1. + (k0 - dk)/2.0) + (ef - 2.0 * Af) * cos_th + (1.0 - 3.0 * (k0 - dk)/2.0) * cos_th*cos_th )
  
def dif_param_comb(cos_th, xs0, Ae, Af, ef, k0, dk):
  """ Combine the two chiral ones with a trick:
      RL cos_th values are shifted by +3 
        => Recognize RL by x>1
  """
  if cos_th>1:
    # Here need to substract the three again!
    return dif_param_RL(cos_th-3, xs0, Ae, Af, ef, k0, dk)
  else:
    return dif_param_LR(cos_th, xs0, Ae, Af, ef, k0, dk)


log.basicConfig(level=log.INFO) # Set logging level
PDF.set_default_mpl_format()
MCLumi = 5000 # MC Statistics is 5ab^-1

input_dir = "/home/jakob/DESY/MountPoints/DUST/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/PrEWInput/MuAcc_costheta_0.9925"

output_dir = input_dir + "/shape_checks"
IOSH.create_dir(output_dir)

LR_RL_pairs = [
  ["2f_mu_81to101_BZ_250_eLpR.csv","2f_mu_81to101_BZ_250_eRpL.csv"],
  ["2f_mu_81to101_FZ_250_eLpR.csv","2f_mu_81to101_FZ_250_eRpL.csv"],
  ["2f_mu_180to275_250_eLpR.csv","2f_mu_180to275_250_eRpL.csv"] ]

for LR_file, RL_file in LR_RL_pairs:
  LR_file_path = input_dir + "/" + LR_file
  RL_file_path = input_dir + "/" + RL_file
  
  log.info("LR file: {}.csv".format(LR_file))
  log.info("RL file: {}.csv".format(RL_file))

  # Read the input files
  LR_reader = IOR.Reader(LR_file_path)
  RL_reader = IOR.Reader(RL_file_path)

  # Get the pandas dataframe for the cut histograms
  angle = "costh_f_star"
  
  LR_df = LR_reader["Data"]
  LR_bin_vals = np.array(LR_df["Cross sections"])
  LR_bin_middles = np.array(LR_df["BinCenters:{}".format(angle)])
  LR_edges_min = np.array(LR_df["BinLow:{}".format(angle)])
  LR_edges_max = np.array(LR_df["BinUp:{}".format(angle)])
  LR_bin_vals *= MCLumi
  
  RL_df = RL_reader["Data"]
  RL_bin_vals = np.array(RL_df["Cross sections"])
  RL_bin_middles = np.array(RL_df["BinCenters:{}".format(angle)])
  RL_edges_min = np.array(RL_df["BinLow:{}".format(angle)])
  RL_edges_max = np.array(RL_df["BinUp:{}".format(angle)])
  RL_bin_vals *= MCLumi
  
  # Shift RL cos_th values for common function trick
  RL_bin_vals += 3
  RL_bin_middles += 3
  RL_edges_min += 3
  RL_edges_max += 3
  
  comb_edges_min = np.concatenate((LR_edges_min,RL_edges_min))
  comb_edges_max = np.concatenate((LR_edges_max,RL_edges_max))
  comb_bin_vals = np.concatenate((LR_bin_vals,RL_bin_vals))

  # Perform fits to the distributions
  fit_vals, p, cov = SST.fit_1D(dif_param_comb, comb_bin_vals, comb_edges_min, comb_edges_max, bounds=[[0,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]])
  
  log.info("Parameters:\n{}".format(p))
  log.info("Uncertainties:\n{}".format(np.sqrt(cov.diagonal())))
  
  # Parameters: xs0, Ae, Af, ef, k0, dk
  i_Ae = 1
  i_Af = 2
  i_ef = 3
  Ae = p[i_Ae]
  Af = p[i_Af]
  ef = p[i_ef]
  AFB = 3/8 * (ef + 2*Ae*Af)
  d_AFB = 3/8 * np.sqrt( \
    cov[i_ef][i_ef] + (2*Ae)**2 * cov[i_Af][i_Af] + (2*Af)**2 * cov[i_Ae][i_Ae] + \
    4 * Ae * cov[i_ef][i_Af] + 4 * Af * cov[i_ef][i_Ae] + 8 * Ae * Af * cov[i_Ae][i_Af]
    )
    
  log.info("AFB: {} \pm {}".format(AFB,d_AFB))
  
  log.info("\n")