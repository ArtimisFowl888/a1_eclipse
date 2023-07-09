import numpy as np


def check_hrrd_gap(zraw, uraw, m_z, m_u, z0, z1, z0p, z1p, note, madz=1000, bndy=0):
    """
  Check the gaps in HRRD data.

  Args:
    zraw: The raw z data.
    uraw: The raw u data.
    m_z: The missing flag for z.
    m_u: The missing flag for u.
    z0: The bottom z (m).
    z1: The top z (m).
    z0p: The extended z0 (m).
    z1p: The extended z1 (m).
    note: The descriptive note on the file being analyzed.
    madz: The maximum allowed value of gap (m). Defaults to 1 km.
    bndy: The maximum gap in the boundary (m). Defaults to 0 km.

  Returns:
    flag: 0 if there are enough data, 1 otherwise.
    z_out: The selected z data.
    u_out: The selected u data.
  """

    flag = 0  # 0: enough data, 1: not enough data
    if not madz:
        madz = 1000
    if not bndy:
        bndy = 0

    z = zraw
    u = uraw

    # (1) pick up z,u that correspond to where z and u are not missing
    ava = np.where(z != m_z and u != m_u)[0]
    if len(ava) == 0:
        print(note + ' No data at all')
        flag = 1
        return
    z = z[ava]
    u = u[ava]

    # (2) pick up z,u where z is within the extended range [z0p,z1p]
    ava = np.where(z >= z0p and z <= z1p)[0]
    if len(ava) == 0:
        print(note + ' No data within the extended altitude range')
        flag = 1
        return
    z = z[ava]
    u = u[ava]

    # (3) check the boundaries, if ztop <= z1-1km or zbot >= z0+1km, return
    if np.max(z) <= z1 - bndy or np.min(z) >= z0 + bndy:
        print(note + ' Not enough boundary coverage')
        flag = 1
        return

    # (4) sort z so that z is in ascending order, excluding same level data
    ava = np.unique(z, return_index=True)[1]
    nz = len(ava)
    nz0 = len(z)
    # if nz != nz0:
    #     print(note+' has '+str(nz0-nz)+' same-z data removed')

    z = z[ava]
    u = u[ava]

    # (5) check the data gap inside, if larger than 1 km, return
    ddzz = z[1:nz - 1] - z[0:nz - 2]
    if np.max(ddzz) >= madz:
        print(note + ' Too large data gap inside')
        flag = 1
        return

    u_out = u
    z_out = z
    return flag, z_out, u_out
