import numpy as np

def phi2kl(phi, l_h):

  """
  Calculate the zonal and meridional wavenumbers from the horizontal phase
  propagation direction and horizontal wavelength.

  Args:
    phi: Horizontal phase propagation direction in degrees.
    l_h: Horizontal wavelength in kilometers.

  Returns:
    k: Zonal wavenumber.
    l: Meridional wavenumber.

  """

  kl = 2.e-3 * np.pi / l_h
  phi2 = phi * np.pi / 180.
  k = kl * np.cos(phi2)
  l = kl * np.sin(phi2)

  return k, l

