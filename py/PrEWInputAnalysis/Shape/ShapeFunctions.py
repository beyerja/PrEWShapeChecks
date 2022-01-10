# ------------------------------------------------------------------------------

""" Functions defining specific shapes that can be tested to fit on a 
    distribution.
"""

# ------------------------------------------------------------------------------

def helicity_amplitudes(x, A_ss, A_os):
  """ Helicity amplitude shape for difermion scattering.
      x ... cos(theta) in difermion rest frame
      A_ss ... same-sign amplitude
      A_os ... opposite-sign amplitude
  """
  return A_ss * (1+x)**2 + A_os * (1-x)**2

def helicity_amplitudes_cor(x, A_ss, A_os, K):
  """ Helicity amplitude shape for difermion scattering including correction 
      term (e.g. for radiative corrections) that leaves integral intact and 
      changes relative ratio of linear to quadratic.
      x ... cos(theta) in difermion rest frame
      A_ss ... same-sign amplitude
      A_os ... opposite-sign amplitude
      K ... correction term factor
  """
  return A_ss * (1+x)**2 + A_os * (1-x)**2 + K * (1 - 3*x**2)