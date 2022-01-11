""" Helper functions related to file name conventions.
"""

import re

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