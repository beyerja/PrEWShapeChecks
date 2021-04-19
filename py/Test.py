import logging as log
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

sys.path.append("./FuncHelp")
import Wrappers as FHW
sys.path.append("./IO")
import Reader as IOR
import SysHelpers as IOSH
sys.path.append("./Shape")
import ShapeFunctions as SSF
import ShapeTesting as SST

log.basicConfig(level=log.INFO) # Set logging level
MCLumi = 5000 # MC Statistics is 5ab^-1

input_dir = "/home/jakob/DESY/MountPoints/DUST/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/PrEWInput"

output_dir = input_dir + "/shape_checks"
IOSH.create_dir(output_dir)

log.info("Looking in dir: {}".format(input_dir))
for file_path in IOSH.find_files(input_dir, ".csv"):
  # Read the input file
  base_name = os.path.basename(file_path).replace(".csv","")
  log.info("Reading file: {}.csv".format(base_name))
  reader = IOR.Reader(file_path)

  # Get the pandas dataframe for the cut histograms
  df = reader["Data"]
  bin_vals = np.array(df["Cross sections"])
  bin_middles = np.array(df["BinCenters:costh_f_star"])
  edges_min = np.array(df["BinLow:costh_f_star"])
  edges_max = np.array(df["BinUp:costh_f_star"])

  # Rescale to MC Lumi
  bin_vals *= MCLumi

  # Perform fits to the distributions
  fit_vals_ha, p_ha, cov_ha = SST.fit_1D(SSF.helicity_amplitudes, bin_vals, edges_min, edges_max, bounds=[[0,0],[np.inf,np.inf]])
  fit_vals_hac, p_hac, cov_hac = SST.fit_1D(SSF.helicity_amplitudes_cor, bin_vals, edges_min, edges_max, bounds=[[0,0,-np.inf],[np.inf,np.inf,np.inf]])
  
  # Calculate the resulting chi-squared's (normalised by degrees of freedom)
  n_bins = len(bin_vals)
  chisq_ndf_ha = SST.chi_squared(bin_vals, fit_vals_ha) / (n_bins - 2)
  chisq_ndf_hac = SST.chi_squared(bin_vals, fit_vals_hac) / (n_bins - 3)
  
  # Before plotting, scale back to cross section level
  bin_vals /= MCLumi
  fit_vals_ha /= MCLumi
  fit_vals_hac /= MCLumi
  yerr = np.sqrt(bin_vals/MCLumi)
  
  # Plot the distributions and the fit results
  fig, ax = plt.subplots(tight_layout=True, figsize=(5,4))
  ax.set_title(base_name)
  ax.errorbar(bin_middles,bin_vals,yerr=yerr, fmt='.',ms=8,mew=2, label="MC")
  ax.plot(bin_middles,fit_vals_ha, label="Pure Hel. Ampl.")
  ax.plot(bin_middles,fit_vals_hac, label="+ corr. term")
  
  ax.set_xlabel(r"$\cos(\theta)$")
  ax.set_ylabel(r"$d\sigma$")
  
  ax.set_ylim(0,ax.get_ylim()[1])

  ax.legend( loc=0, title=r"$\chi^2_{pure}/ndf = $" + str(np.round_(chisq_ndf_ha, decimals=1)) + "\n$\chi^2_{cor}/ndf = $" + str(np.round_(chisq_ndf_hac, decimals=1)) )

  fig.savefig("{}/{}_shape_check.pdf".format(output_dir, base_name))
