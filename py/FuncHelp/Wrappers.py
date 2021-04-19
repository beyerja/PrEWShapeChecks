# ------------------------------------------------------------------------------

""" Useful function wrappers.
"""

# ------------------------------------------------------------------------------

def array_arg_wrapper(func, args):
  """ Allow calling a function with an array of its arguments.
  """
  return func(*args)
