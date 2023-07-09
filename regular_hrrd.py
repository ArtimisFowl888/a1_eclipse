import numpy as np


def regular_hrrd(zbig, z, u, flag):
    """
  Interpolate HRRD data to regular grid.

  Args:
    zbig: The extended z-grid (m).
    z: The raw z data.
    u: The raw u data.
    flag: Whether there are enough valid data.

  Returns:
    ugrid: The interpolated u data on the regular grid.
    zgrid: The corresponding regular z-grid.
  """

    nodata = -999.

    if flag == 1:
        ugrid = nodata
        zgrid = nodata
        return

    # Find the indices of zgrid that correspond to z
    i0 = np.where(zbig >= min(z))[0][0]
    i1 = np.where(zbig >= max(z))[0][0] - 1

    # Interpolate u to the regular grid
    zgrid = zbig[i0:i1]
    ugrid = np.interp(zgrid, z, u)

    return ugrid, zgrid
