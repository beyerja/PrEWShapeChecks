# ------------------------------------------------------------------------------

""" Helper function related to the system that your running on.
"""

# ------------------------------------------------------------------------------

from pathlib import Path
import os

# ------------------------------------------------------------------------------
# Creating

def create_dir(dir):
  """ Try to create the given directory and it's parents.
  """
  Path(dir).mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------------------
# Finding

def find_files(dir,extension):
  """ Find all the files of a given extension in the directory.
  """
  return["{}/{}".format(dir,dir_object) for dir_object in os.listdir(dir) if extension in dir_object]

# ------------------------------------------------------------------------------
