import pandas as pd

# Local modules
import IO.CSVMetadataReader as CMR

# ------------------------------------------------------------------------------

class Reader:
  """ Class to read / interpret a distribution data file.
  """
  
  # --- Constructor ------------------------------------------------------------
  
  def __init__(self,file_path):
    self.data = {}
    self.interpret(file_path)
    
  # --- Access functions -------------------------------------------------------
    
  def __getitem__(self,index):
    """ Define what happens when the [] operator is applied (only reading).
    """
    return self.data[index]
    
  # --- Internal functions -----------------------------------------------------
    
  def interpret(self,file_path):
    """ Read and interpret the input file.
    """
    # Find and use the metadata
    mr = CMR.CSVMetadataReader(file_path)
    self.data = mr.metadata
    
    # Find the distribution data
    self.data["Data"] = pd.read_csv(file_path, header=mr.get_data_header_line())
    
# ------------------------------------------------------------------------------