import numpy as np

def get_coriolis(latitude):
  """
  Get the Coriolis parameter at a given latitude.

  Args:
    latitude: The latitude in degrees.

  Returns:
    The Coriolis parameter in s^-1.
  """

  omega = 7.2919996e-05
  f = 2.0 * omega * np.sin(latitude * np.pi / 180.)
  return f

