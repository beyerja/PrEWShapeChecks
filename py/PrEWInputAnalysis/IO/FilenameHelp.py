""" Helper functions related to file name conventions.
"""

import re

def find_Z_direction(file_path):
  """ Find the Z direction in the file path.
  """
  if "_BZ_" in file_path:
    return "BZ"
  elif "_FZ_" in file_path:
    return "FZ"
  else:
    return None
    
def find_2f_mass_label(file_path):
  """ Find the mass label of the 2f process in the file path.
  """
  if "180to275" in file_path:
    return "high_Q2"
  elif "81to101" in file_path:
    return "return_to_Z"
  else:
    raise Exception("Could not find 2f mass label in ", file_path)

def find_chirality(file_path):
  """ Find the chirality in the given file path.
  """
  eM_chi = re.search(r"e(L|R)",file_path).group(0)
  eP_chi = re.search(r"p(L|R)",file_path).group(0)
  return "{}{}".format(eM_chi,eP_chi)

def find_lep_charge(file_path, lep):
  """ Find the chirality in the given file path.
  """
  if lep+"minus" in file_path:
    return -1
  elif lep+"plus" in file_path:
    return +1
  else:
    raise Exception("Did not find {} charge in {}".format(lep, file_path))