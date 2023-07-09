import numpy as np


def perturb_hrrd(u, flag, npoly=None, flow=None, fhigh=None):
    """
  Perturb HRRD data.

  Args:
    u: The interpolated u data on the regular grid.
    flag: Whether there are enough valid data.
    npoly: The order of polynomial fit. Defaults to None.
    flow: The lower limit of the high-pass filter. Defaults to None.
    fhigh: The upper limit of the low-pass filter. Defaults to None.

  Returns:
    up: The perturbed u data.
    ub: The background u data.
  """

    if flag == 1:
        return

    # Polynomial fit
    if npoly is not None:
        co = np.polyfit(np.arange(len(u)), u, npoly, full=True)
        up = u - co[0]
        return up, co[0]

    # High-pass filter
    if flow is not None:
        up = np.digital_smooth(u, flow, 1)
        ub = u - up
        return up, ub

    # Low-pass filter
    if fhigh is not None:
        ub = np.digital_smooth(u, 0, fhigh)
        up = u - ub
    return up, ub
