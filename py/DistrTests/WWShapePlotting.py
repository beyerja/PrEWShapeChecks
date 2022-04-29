import logging as log
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from tqdm import tqdm

sys.path.append("../PrEWInputAnalysis")
import FuncHelp.Wrappers as FHW
import IO.FilenameHelp as IOFH
import IO.Reader as IOR
import IO.SysHelpers as IOSH
import Plotting.DefaultFormat as PDF
import Plotting.Naming as PN

def add_hist2d(ax, x, y, nbins, xmin, xmax, i, j, **kwargs):
  """ Add the hist2d to the ax for the angles of indices i and j.
  """
  return ax.hist2d(
    x=x[i], y=x[j], weights=y,
    bins=(nbins[i], nbins[j]), range=[(xmin[i], xmax[i]),(xmin[j], xmax[j])],
    **kwargs)

def create_2D_projection_plot(x, y, xmin, xmax, nbins, angles, base_name, 
                              outdir, out_formats):
  """ Create a plot showing the 3 different 2D projections of the 3D 
      distribution.
  """
  # Figure basics
  fig, axs = plt.subplots(2, 3, sharex='col', sharey='row', figsize=(12,9),
                          tight_layout=True,
                          gridspec_kw={'hspace': 0.0, 'wspace': 0.0})
  (ax1, ax_empty, _1), (ax2, ax3, _2) = axs
  _1.axis('off')
  _2.axis('off')
  ax_empty.axis('off')

  # Draw ticks only on the relevant axes
  ax1.tick_params(bottom=False, top=True, left=True, right=True, 
                  labelleft=False, labeltop=True, labelright=True)
  ax2.tick_params(labelbottom=False, labelleft=False)
  ax3.tick_params(bottom=True, top=True, left=False, right=True, 
                  labeltop=True, labelright=True, labelbottom=False)

  # Actually plot the histograms
  h1 = add_hist2d(ax1, x, y, nbins, xmin, xmax, 0, 1)
  h2 = add_hist2d(ax2, x, y, nbins, xmin, xmax, 0, 2)
  h3 = add_hist2d(ax3, x, y, nbins, xmin, xmax, 1, 2)

  # Set a useful common color scale on all histograms
  cmax = np.amax([ np.amax(y_h) for y_h in (h1[0], h2[0], h3[0]) ])
  for h in [h1,h2,h3]:
    h[-1].set_clim(0, cmax)

  # Set useful axis limits
  ax1.set_xlim((xmin[0],xmax[0]))
  ax1.set_ylim((xmin[1],xmax[1]))
  ax3.set_xlim((xmin[1],xmax[1]))
  ax3.set_ylim((xmin[2],xmax[2]))

  # Create small white lines between the plots
  for ax_line in (ax2.spines['top'],ax2.spines['right'],ax1.spines['bottom'],
                  ax3.spines['left']):
    ax_line.set_linewidth(3)
    ax_line.set_color('white')

  # Proper labelling of the axes
  label_args = {"fontsize":30, "ha":'center'}
  ax1.text(0.5, 1.2, "${}$".format(PN.observable_str(angles[0], "WW")), transform=ax1.transAxes, **label_args)
  ax_empty.text(0.5, 0.5, "${}$".format(PN.observable_str(angles[1], "WW")), transform=ax_empty.transAxes, **label_args)
  ax3.text(1.12, 0.65, "${}$".format(PN.observable_str(angles[2], "WW")), transform=ax3.transAxes, **label_args)

  # Add a colorbar
  # cbar_title = "$\\frac{d^2\\sigma}{dx_1dx_2} [$fb$^{{-1}}]$"
  cbar_title = "$d\\sigma [$fb$]$"
  ax_whole = fig.add_subplot(1,6,5, visible=False)
  fig.colorbar(h1[-1], ax=ax_whole, label=cbar_title, fraction=1.0)
  # ax_whole.text(0.5, 0.65, "$d\\sigma [$fb$^{{-1}}]$", 
  #               transform=ax_whole.transAxes, **label_args)

  # Mark which process it is
  chirality = IOFH.find_chirality(base_name)
  mu_charge = IOFH.find_lep_charge(base_name, "mu")
  process_str="${}$".format(PN.WW_process_str(chirality, mu_charge))
  ax_empty.text(0.5, 1.2, process_str, transform=ax_empty.transAxes, 
        bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'), 
        **label_args)

  # Save the plot in files
  for out_format in out_formats:
    format_dir = "{}/{}".format(outdir,out_format)
    IOSH.create_dir(format_dir)
    fig.savefig("{}/{}.{}".format(format_dir,base_name,out_format), 
                transparent=True, bbox_inches='tight')
  plt.close(fig)
  
def create_1D_projection_plot(i_coord, x, y, xmin, xmax, nbins, angles, 
                              base_name, outdir, out_formats):
  """ Create a 1D projection plot for the given coordinate.
  """
  
  # Create the figure and plot
  fig = plt.figure(figsize=(6.5, 5), tight_layout=True)
  ax = plt.gca()
  ax.hist(
    x=x[i_coord], weights=y,
    bins=nbins[i_coord], range=(xmin[i_coord], xmax[i_coord]),
    ls="-", lw=3, histtype=u'step')
    
  # Useful limits
  ax.set_xlim(xmin[i_coord], xmax[i_coord])
  ax.set_ylim(0, 1.1*ax.get_ylim()[1])
  
  # Useful labels
  angle = angles[i_coord]
  ax.set_xlabel("${}$".format(PN.observable_str(angle, "WW")))
  ax.set_ylabel("$d\\sigma [$fb$]$")
  
  # Mark which process it is
  chirality = IOFH.find_chirality(base_name)
  mu_charge = IOFH.find_lep_charge(base_name, "mu")
  process_str="${}$".format(PN.WW_process_str(chirality, mu_charge))
  ax.set_title(process_str)
  
  # Save the plot in files
  for out_format in out_formats:
    format_dir = "{}/{}".format(outdir,out_format)
    IOSH.create_dir(format_dir)
    fig.savefig("{}/{}_{}.{}".format(format_dir,base_name,angle,out_format), 
                transparent=True, bbox_inches='tight')

def plot_WW_distr(infile, outdir, out_formats=["pdf","png"]):
  """ Plot the distribution that is stored in the given file.
  """
  base_name = os.path.basename(infile).replace(".csv","")
  
  log.debug("Reading file: {}.csv".format(base_name))
  reader = IOR.Reader(infile)
  
  # Get the pandas dataframe for the cut histograms
  df = reader["Data"]
  
  # Get the data
  angles = ("costh_Wminus_star", "costh_l_star", "phi_l_star")
  x = np.array([df["BinCenters:{}".format(angle)] for angle in angles])
  y = np.array(df["Cross sections"])
  
  xmin = np.array([np.amin(df["BinLow:{}".format(angle)]) for angle in angles])
  xmax = np.array([np.amax(df["BinUp:{}".format(angle)]) for angle in angles])
  
  # Bin counting assumes that bins are not binned thinner than 1/1000 and that
  # absolute values are of order 0.1-1
  nbins = np.array([
    len( np.unique((x[d]*10000.).astype(int)) ) for d in range(len(x)) ])
  
  create_2D_projection_plot(
    x, y, xmin, xmax, nbins, angles, base_name, outdir, out_formats)
  
  for i in range(len(angles)):
    create_1D_projection_plot(i, x, y, xmin, xmax, nbins, angles, base_name, 
                              outdir, out_formats)
  
def main():
  log.basicConfig(level=log.INFO) # Set logging level
  PDF.set_default_mpl_format()
  input_dir = "/home/jakob/DESY/MountPoints/DUST/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/PrEWInput"
  output_dir = input_dir + "/shape_checks/WWShapePlots"
  
  log.info("Looking in dir: {}".format(input_dir))
  for file_path in tqdm(IOSH.find_files(input_dir, ".csv"), desc="files"):
    if "tau" in file_path:
      continue # Skipping tau distributions, not use right now
    plot_WW_distr(file_path, output_dir)

if __name__ == "__main__":
  main()