import numpy as np


def dispersion_m_w_to_k(w, m, lat, n=0.02, hrou=7., help=False):
  """
  Calculate the horizontal wavenumber from the intrinsic frequency and vertical
  wavenumber.

  Args:
    w: The intrinsic frequency in s^-1.
    m: The vertical wavenumber in s^-1.
    lat: The latitude in degrees.
    n: The Brunt-Väisälä frequency in s^-1.
    hrou: The density scale height in km.
    help: If True, print help information.

  Returns:
    The horizontal wavenumber in s^-1.
  """

  if help:
    print("Purpose:")
    print("   ( intrinsic frequency + vertical wavenumber)")
    print("                 --> horizontal wavenumber")
    print("Usage:")
    print("    k = dispersion_m_w_to_k(w, m, f, n=n, alpha=alpha)")
    print("Steps:")
    print("   ")
    print("Input:")
    print("    w -- intrinsic frequency [SI]")
    print("    m -- vertical wavenumber [SI], 2*!PI / Lz")
    print("    lat -- latitude [deg]")
    print("Keywords:")
    print("    help -- print this information")
    print("    n -- Brunt-Vasalla freq [SI]")
    print("    hrou -- density scale height [km]")
    print("Output:")
    print("    k -- horizontal wavenumber [SI], 2*!PI / Lh")
    print("Examples:")
    print("   ")
    print("History:")
    print("    08/11/2006 created by Jie Gong")
    print("Note:")
    print("   ")
    print("")
    raise SystemExit

  if hrou is None:
    hrou = 7.

  if n is None:
    n = 0.02  # typical stratospheric value

  alpha = 0.5 / (hrou * 1e3)
  f = get_coriolis(lat)

  k = np.sqrt((w * w - f * f) / (n * n - w * w) * (m * m + alpha * alpha))

  return k
