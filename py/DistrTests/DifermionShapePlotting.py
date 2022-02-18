import logging as log
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from tqdm import tqdm

sys.path.append("../PrEWInputAnalysis")
import IO.FilenameHelp as IOFH
import IO.Reader as IOR
import IO.SysHelpers as IOSH
import Plotting.DefaultFormat as PDF
import Plotting.Naming as PN

def plot_mumu_distr(infile, outdir, out_formats=["pdf","png"]):
  """ Plot the distribution that is stored in the given file.
  """
  base_name = os.path.basename(infile).replace(".csv","")
  
  log.debug("Reading file: {}.csv".format(base_name))
  reader = IOR.Reader(infile)
  
  # Get the pandas dataframe for the cut histograms
  df = reader["Data"]
  
  # Get the data
  x_name = "costh_f_star"
  x = np.array(df["BinCenters:{}".format(x_name)])
  y = np.array(df["Cross sections"])
  
  xmin = np.array(np.amin(df["BinLow:{}".format(x_name)]))
  xmax = np.array(np.amax(df["BinUp:{}".format(x_name)]))
  
  # Bin counting assumes that bins are not binned thinner than 1/1000 and that
  # absolute values are of order 0.1-1
  nbins = len( np.unique((x*10000.).astype(int)) )
  
  # Create the figure and plot
  fig = plt.figure(figsize=(6.5, 5), tight_layout=True)
  ax = plt.gca()
  ax.hist( x=x, weights=y, bins=nbins, range=(xmin, xmax),
           ls="-", lw=3, histtype=u'step' )
    
  # Useful limits
  ax.set_xlim(xmin, xmax)
  # ax.set_ylim(0, 1.1*ax.get_ylim()[1])
  
  # Useful labels
  ax.set_xlabel("${}$".format(PN.observable_str(x_name, "mumu")))
  ax.set_ylabel("$d\\sigma [$fb$^{{-1}}]$")
  
  # Mark which process it is
  chirality = IOFH.find_chirality(base_name)
  Z_direction = IOFH.find_Z_direction(base_name)
  mass_label = IOFH.find_2f_mass_label(base_name)
  process_str = "${}$".format(PN.difermion_process_str(
                               "mu", chirality, mass_label, Z_direction))
  ax.set_title(process_str)
  
  # Save the plot in files
  for out_format in out_formats:
    format_dir = "{}/{}".format(outdir,out_format)
    IOSH.create_dir(format_dir)
    fig.savefig("{}/{}_{}.{}".format(format_dir,base_name,x_name,out_format), 
                transparent=True, bbox_inches='tight')
                
def main():
  log.basicConfig(level=log.INFO) # Set logging level
  PDF.set_default_mpl_format()
  input_dir = "/home/jakob/DESY/MountPoints/DUST/TGCAnalysis/SampleProduction"+\
              "/NewMCProduction/2f_Z_l/PrEWInput/MuAcc_costheta_0.9925"
  output_dir = input_dir + "/shape_checks/2fShapePlots"
  
  log.info("Looking in dir: {}".format(input_dir))
  for file_path in tqdm(IOSH.find_files(input_dir, ".csv"), desc="files"):
    if not "2f_mu" in file_path:
      continue # Only draw mumu distributions
    plot_mumu_distr(file_path, output_dir)

if __name__ == "__main__":
  main()