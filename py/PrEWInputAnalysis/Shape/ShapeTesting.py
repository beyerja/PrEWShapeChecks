# ------------------------------------------------------------------------------

""" Functions and classes to test if a given distributions follows a specified 
    shape.
"""

# ------------------------------------------------------------------------------

# External packages
import functools
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import sys

# Local packages
import FuncHelp.Wrappers as FHW

# ------------------------------------------------------------------------------

def bin_integral_1D(func):
  @functools.wraps(func)
  def wrapper_bin_integrated(x, *args):
    """ Return the integral of the given function in each bin.
    """
    f = lambda _x: func(_x, *args)
    return np.array([integrate.quad(f, x[b][0], x[b][1])[0] for b in range(len(x))])
  return wrapper_bin_integrated

# ------------------------------------------------------------------------------

def fit_1D(func, bin_vals, edges_min, edges_max, p0=None, bounds=(-np.inf,np.inf)):
  """ Fit the given function to the provided values (which are the integral in 
      each bin).
  """
  x = np.column_stack((edges_min, edges_max)) # x ... Edges
  yerr = np.sqrt(bin_vals) # Gaussian error
  f = bin_integral_1D(func)
  p, cov = curve_fit(f=f, xdata=x, ydata=bin_vals, sigma=yerr, p0=p0, bounds=bounds)
  fit_y = FHW.array_arg_wrapper(f,[x, *p]) 
  
  return fit_y, p, cov

# ------------------------------------------------------------------------------

def chi_squared(bin_vals, fit_vals):
  return np.sum( (bin_vals - fit_vals)**2 / bin_vals )

# ------------------------------------------------------------------------------

