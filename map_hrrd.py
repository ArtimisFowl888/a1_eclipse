import numpy as np

def map_hrrd(flag, zgrid, zm, up, ub):
  """
  Map HRRD data to a common z-grid

  Args:
    flag: 0: enough valid data; 1: not enough
    zgrid: extended z-grid (m) specific to the given variable
    zm: common z-grid (desired)
    up: perturbation
    ub: background

  Returns:
    up: perturbation on common z-grid
    ub: background on common z-grid
  """

  nz = len(zm)

  # Not enough data

  if flag == 1:
    up = np.full(nz, -999.)
    ub = np.full(nz, -999.)
    return up, ub

  # Find the indices of the z-grid points that are within the range of zm

  i0 = np.where(zgrid >= min(zm))[0][0]
  i1 = np.where(zgrid >= max(zm))[0][0] - 1

  # Return the perturbation and background on the common z-grid

  up = up[i0:i1]
  ub = ub[i0:i1]

  return up, ub

