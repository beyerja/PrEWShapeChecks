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

input_dir_nocor = "/home/jakob/DESY/MountPoints/DUST/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/PrEWInput/MuAcc_costheta_0.9925"
input_dir_sicor = "/home/jakob/DESY/MountPoints/DUST/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/PrEWInput/MuAcc_costheta_0.9925/TrueAngle"

output_dir = input_dir_nocor + "/ISR_shape_checks"
IOSH.create_dir(output_dir)

file_pairs = [
  ["2f_mu_81to101_BZ_250_eLpR.csv", "2f_mu_81to101_BZ_true_250_eLpR.csv"],
  ["2f_mu_81to101_BZ_250_eRpL.csv", "2f_mu_81to101_BZ_true_250_eRpL.csv"],
  ["2f_mu_81to101_FZ_250_eLpR.csv", "2f_mu_81to101_FZ_true_250_eLpR.csv"],
  ["2f_mu_81to101_FZ_250_eRpL.csv", "2f_mu_81to101_FZ_true_250_eRpL.csv"],
]

for file_pair in file_pairs:
  # Read the input file
  base_name = file_pair[0].replace(".csv","")
  log.info("Processing: {}".format(base_name))
  reader_nocor = IOR.Reader(input_dir_nocor + "/" + file_pair[0])
  reader_sicor = IOR.Reader(input_dir_sicor + "/" + file_pair[1])

  # Get the pandas dataframe for the cut histograms
  df_nocor = reader_nocor["Data"]
  df_sicor = reader_sicor["Data"]
  bin_vals_nocor = MCLumi * np.array(df_nocor["Cross sections"])
  bin_vals_sicor = MCLumi * np.array(df_sicor["Cross sections"])
  
  angle = "costh_f_star"
  bin_middles = np.array(df_nocor["BinCenters:{}".format(angle)])
  edges_min = np.array(df_nocor["BinLow:{}".format(angle)])
  edges_max = np.array(df_nocor["BinUp:{}".format(angle)])
  bin_width = edges_max[0] - edges_min[0]

  # Perform fits to the distributions
  fit_vals_nocor, p_nocor, cov_nocor = SST.fit_1D(SSF.helicity_amplitudes, bin_vals_nocor, edges_min, edges_max, bounds=[[0,0],[np.inf,np.inf]])
  fit_vals_sicor, p_sicor, cov_sicor = SST.fit_1D(SSF.helicity_amplitudes, bin_vals_sicor, edges_min, edges_max, bounds=[[0,0],[np.inf,np.inf]])
  
  # log.info("Parameters:")
  # log.info("  Pure HA: {}".format(p_ha))
  # log.info("  w/ corr: {}".format(p_hac))
  
  # Calculate the resulting chi-squared's (normalised by degrees of freedom)
  n_bins = len(bin_vals_nocor)
  chisq_ndf_nocor = SST.chi_squared(bin_vals_nocor, fit_vals_nocor) / (n_bins - 2)
  chisq_ndf_sicor = SST.chi_squared(bin_vals_sicor, fit_vals_sicor) / (n_bins - 3)
  
  # Before plotting, scale back to cross section level
  bin_vals_nocor /= MCLumi * bin_width
  bin_vals_sicor /= MCLumi * bin_width
  fit_vals_nocor /= MCLumi * bin_width
  fit_vals_sicor /= MCLumi * bin_width
  yerr_nocor = np.sqrt(bin_vals_nocor/MCLumi)
  yerr_sicor = np.sqrt(bin_vals_sicor/MCLumi)
  
  # --- Plot both distributions ------------------------------------------------
  
  # Plot the distributions and the fit results
  fig, ax = plt.subplots(tight_layout=True, figsize=(8.5,7))
  ax.plot(bin_middles,fit_vals_nocor, ls="--", lw=2, label="Fit no cor.", zorder=1, alpha=0.7)
  ax.plot(bin_middles,fit_vals_sicor, ls="--", lw=2, label="Fit w/ cor.", zorder=1, alpha=0.7)
  ax.errorbar(bin_middles,bin_vals_nocor,yerr=yerr_nocor, fmt='.',ms=12,mew=2, label="no ISR cor.", zorder=2)
  ax.errorbar(bin_middles,bin_vals_sicor,yerr=yerr_sicor, fmt='.',ms=12,mew=2, label="w/ ISR cor.", zorder=2)
  
  ax.set_xlabel("${}$".format(PN.observable_str(angle, "mumu")))
  ax.set_ylabel(r"$\frac{d\sigma}{d\,\cos\theta} [$fb$]$")
  
  ax.set_xlim(-1,1)
  ax.set_ylim(0,0.98*ax.get_ylim()[1])
  
  # Mark which process it is
  chirality = IOFH.find_chirality(base_name)
  Z_direction = IOFH.find_Z_direction(base_name)
  mass_label = IOFH.find_2f_mass_label(base_name)
  process_str = "${}$".format(PN.difermion_process_str(
                               "mu", chirality, mass_label, Z_direction))
  ax.set_title(process_str)

  handles, _ = ax.get_legend_handles_labels()
  reordering = [2,3,0,1]
  ax.legend( handles=[handles[i] for i in reordering], loc=0, ncol=2, 
             title=r"$\chi^2_{no\, cor.}\,/ndf = $" + \
                    str(np.round_(chisq_ndf_nocor, decimals=1)) + \
                    "\n$\chi^2_{w/\, cor.}\,/ndf = $" + \
                    str(np.round_(chisq_ndf_sicor, decimals=1)) )

  fig.savefig("{}/{}_ISR_shape_check.pdf".format(output_dir, base_name))
  
  plt.close(fig)
