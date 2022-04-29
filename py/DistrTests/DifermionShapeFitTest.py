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

log.basicConfig(level=log.INFO) # Set logging level
PDF.set_default_mpl_format()
MCLumi = 5000 # MC Statistics is 5ab^-1

input_dir = "/home/jakob/DESY/MountPoints/DUST/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/PrEWInput/MuAcc_costheta_0.9925"
# input_dir = "/home/jakob/DESY/MountPoints/DUST/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/PrEWInput/MuAcc_costheta_0.9925/TrueAngle"

output_dir = input_dir + "/shape_checks"
IOSH.create_dir(output_dir)

log.info("Looking in dir: {}".format(input_dir))
for file_path in IOSH.find_files(input_dir, ".csv"):
  # Read the input file
  base_name = os.path.basename(file_path).replace(".csv","")
  log.info("Reading file: {}.csv".format(base_name))
  reader = IOR.Reader(file_path)

  # Get the pandas dataframe for the cut histograms
  angle = "costh_f_star_true" if "TrueAngle" in input_dir else "costh_f_star"
  df = reader["Data"]
  bin_vals = np.array(df["Cross sections"])
  bin_middles = np.array(df["BinCenters:{}".format(angle)])
  edges_min = np.array(df["BinLow:{}".format(angle)])
  edges_max = np.array(df["BinUp:{}".format(angle)])
  bin_width = edges_max[0] - edges_min[0]

  # Rescale to MC Lumi
  bin_vals *= MCLumi

  # Perform fits to the distributions
  fit_vals_ha, p_ha, cov_ha = SST.fit_1D(SSF.helicity_amplitudes, bin_vals, edges_min, edges_max, bounds=[[0,0],[np.inf,np.inf]])
  fit_vals_hac, p_hac, cov_hac = SST.fit_1D(SSF.helicity_amplitudes_cor, bin_vals, edges_min, edges_max, bounds=[[0,0,-np.inf],[np.inf,np.inf,np.inf]])
  
  log.info("Parameters:")
  log.info("  Pure HA: {}".format(p_ha))
  log.info("  w/ corr: {}".format(p_hac))
  
  # Calculate the resulting chi-squared's (normalised by degrees of freedom)
  n_bins = len(bin_vals)
  chisq_ndf_ha = SST.chi_squared(bin_vals, fit_vals_ha) / (n_bins - 2)
  chisq_ndf_hac = SST.chi_squared(bin_vals, fit_vals_hac) / (n_bins - 3)
  
  # Before plotting, scale back to cross section level
  bin_vals /= MCLumi / bin_width
  fit_vals_ha /= MCLumi / bin_width
  fit_vals_hac /= MCLumi / bin_width
  yerr = np.sqrt(bin_vals/MCLumi)
  
  # --- Plot both distributions ------------------------------------------------
  
  # Plot the distributions and the fit results
  fig, ax = plt.subplots(tight_layout=True, figsize=(8.5,7))
  ax.set_title(base_name)
  ax.errorbar(bin_middles,bin_vals,yerr=yerr, fmt='.',ms=12,mew=2, label="MC")
  ax.plot(bin_middles,fit_vals_ha, ls="--", lw=2, label="Pure Hel. Ampl.")
  ax.plot(bin_middles,fit_vals_hac, ls="--", lw=2, label="+ corr. term")
  
  ax.set_xlabel(r"$\cos(\theta)$")
  ax.set_ylabel(r"$\frac{d\sigma}{d\,\cos\theta} [$fb$]$")
  if "TrueAngle" in input_dir:
    ax.set_xlabel(r"$\cos(\theta_{ISR-corrected})$")

  # Mark which process it is
  chirality = IOFH.find_chirality(base_name)
  Z_direction = IOFH.find_Z_direction(base_name)
  mass_label = IOFH.find_2f_mass_label(base_name)
  process_str = "${}$".format(PN.difermion_process_str(
                               "mu", chirality, mass_label, Z_direction))
  ax.set_title(process_str)

  ax.set_xlim(-1,1)
  ax.set_ylim(0,ax.get_ylim()[1])

  ax.legend( loc=0, ncol=2, title=r"$\chi^2_{pure}/ndf = $" + str(np.round_(chisq_ndf_ha, decimals=1)) + "\n$\chi^2_{cor}/ndf = $" + str(np.round_(chisq_ndf_hac, decimals=1)) )

  fig.savefig("{}/{}_shape_check.pdf".format(output_dir, base_name))
  
  plt.close(fig)
  
  # --- Plot only helicity amplitude approach ----------------------------------
  
  # Plot the distributions and the fit results
  fig, ax = plt.subplots(tight_layout=True, figsize=(8.5,7))
  ax.errorbar(bin_middles,bin_vals,yerr=yerr, fmt='.',ms=12,mew=2, label="MC")
  ax.plot(bin_middles,fit_vals_ha, ls="--", lw=2, label="Pure Hel. Ampl.")
  
  ax.set_xlabel(r"$\cos(\theta)$")
  ax.set_ylabel(r"$\frac{d\sigma}{d\,\cos\theta} [$fb$]$")
  
  # Mark which process it is
  chirality = IOFH.find_chirality(base_name)
  ax.set_title(process_str)
  
  ax.set_xlim(-1,1)
  ax.set_ylim(0,ax.get_ylim()[1])

  ax.legend( loc=0, title=r"$\chi^2_{pure}/ndf = $" + str(np.round_(chisq_ndf_ha, decimals=1)) )

  fig.savefig("{}/{}_shape_check_HelAmplOnly.pdf".format(output_dir, base_name))
  plt.close(fig)
