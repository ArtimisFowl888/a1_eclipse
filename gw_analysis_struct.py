import numpy as np

def gw_analysis_struct(input, plot=False, ttl=None):
  """
  Do the GW analysis following VAE97

  Args:
    input: The input structure
    plot: Whether to make a preliminary plot
    ttl: The title of the preliminary plot

  Returns:
    gw: All the GW info. derived from this analysis
  """

  if plot:
    print("Making preliminary plot...")

  ke, kez, kev, kevert, pe, f_w, m_bar, up_frac, mean_dir, uw, vw, n_bar = gw_analysis2(
      input.up, input.vp, input.arp, input.tp, input.um, input.vm, input.arm,
      input.tm, input.pm, input.lat[0], input.lon[0], input.z, input.dz,
      input.ipw, input.iw, input.id)

  lz = 1e-3 / m_bar  # [km]
  m = 2 * np.pi * m_bar  # 2*!PI / Lz instead of 1/Lz
  lat = input.lat[0]  # [deg]
  f = get_coriolis(lat)
  omega = f_w * f  # intrinsic frequency [SI]

  k = dispersion_m_w_to_k(omega, m, lat, n=n_bar, hrou=7.)
  lh = 2e-3 * np.pi / k  # horizontal wavelength [km]
  phi2kl(mean_dir, lh, kx, ky)

  gw = {
      "ke": ke,
      "kez": kez,
      "kev": kev,
      "kevert": kevert,
      "pe": pe,
      "f_w": f_w,
      "m_bar": m_bar,
      "up_frac": up_frac,
      "mean_dir": mean_dir,
      "uw": uw,
      "vw": vw,
      "n_bar": n_bar,
      "lz": lz,
      "lh": lh,
      "m": m,
      "f": f,
      "lat": lat,
      "omega": omega,
      "k": k,
      "kx": kx,
      "ky": ky,
      "input": input,
  }

  return gw

