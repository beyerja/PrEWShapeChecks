import logging as log
import numpy as np

# ------------------------------------------------------------------------------

class CSVMetadataReader:
  """ Class to read the metadata from the top of an CSV file.
  """
  begin_marker = "#BEGIN-METADATA"
  end_marker = "#END-METADATA"
  
  # Arrays storing how specific metadata lines should be interpreted
  # (The ones that aren't mentioned are interpreted as string)
  int_data = ["Energy", "e-Chirality", "e+Chirality"]
  float_data = ["Coef|MuonAcc_CutValue"]

  def __init__(self, csv_path):
    self.metadata = {}
    self.data_header_line = None
    self.interpret(csv_path)

  # --- Access functions -------------------------------------------------------

  def __getitem__(self,index):
    """ Define what happens when the [] operator is applied (only reading).
    """
    return self.metadata[metadata_ID]
      
  def get_data_header_line(self):
    """ Get the line index at which the actual CSV data starts.
    """
    return self.data_header_line

  # --- Internal functions -----------------------------------------------------
    
  def find_metadata_lines(self, csv_path):
    """ Find the lines in the given file that correspond to the metadata.
    """
    # Look line by line
    line_index = 0
    in_metadata = False
    metadata_lines = [] 
    with open(csv_path, 'r') as read_obj:
      # Check line by line
      for line in read_obj:
        line_index += 1
        line = line.strip() # Remove trailing/leading whitespaces etc.
        if self.end_marker in line:
          break
        elif self.begin_marker in line:
          in_metadata = True
        elif in_metadata:
          metadata_lines.append(line)
        else:
          raise ValueError("Unexpected line before CSV Metadata: {}".format(line))
      
    self.data_header_line = line_index
    return metadata_lines
    
  def interpret_metadata(self, ID, data_str):
    """ Interpret the string of a specific metadata.
    """
    data = data_str
    
    if ID in CSVMetadataReader.int_data:
      data = int(float(data))
    elif ID in CSVMetadataReader.float_data:
      data = float(data)
      
    return data
    
  def interpret(self, csv_path):
    """ Interpret the metadata lines found in the CSV file.
    """
    # Split each line by the ":" separator and store its value
    for line in self.find_metadata_lines(csv_path):
      log.debug("Interpreting line: {}".format(line))
      ID, data_str = line.split(":")
      data_str = data_str.strip() # Remove trailing/leading whitespaces etc.
    
      self.metadata[ID] = self.interpret_metadata(ID,data_str)

# ------------------------------------------------------------------------------
