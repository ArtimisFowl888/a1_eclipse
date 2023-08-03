import numpy as np
from scipy import signal

# assumed function imports - the user will need to replace with correct functions
from idl_functions import avg1d, constants, detrend, t_spec_new, t_spec2, postdark, fit, fit2, fit3, fit4, fit5, fit6, fit10, fit11, rotary_spec, directions_new, stokes_new, phase_speed_new, phase_speed2

def gw_analysis2(up, vp, arp, tp, um, vm, arm, tb, pr, lat, lon, zkm, dz, ipw, iw, id, n_bar=None, plot=None, ttl=None):

    fit6,t6_0 = None,None # assuming global variables
    nnn = len(zkm)
    dzm = dz*1000
    zm = zkm*1000
    z = zkm
    p = pr*100
    p_r = p/287.05

    tph = tp/tb
    tph2 = tph.copy()
    dt = np.gradient(tb)/dz
    n_sq = (dt+9.8)*9.8/(tb*1000)
    n_bar = np.sqrt(avg1d(n_sq))

    constants, lat, n_bar, f, bo, co, w1, w2, w1p, deltaw = None,None,None,None,None,None,None,None,None,None
    tps_bar = avg1d(tph**2)

    if id == 1:
        tphp = tph.copy()
        detrend(tphp, z, 1)
        tph2 = tphp

    if ipw == 1:
        tph2 = tph[1:nnn]-tph[0:nnn-1]

    if iw == 0:
        t_spec_new(tph2, tsh, terr, m, lamda, 1, dzm)

    if iw == 1:
        t_spec2(tph2, tsh, terr, m, lamda, 1, dzm)

    if ipw == 1:
        postdark(tsh, m)

    if plot is not None:
        tspec_plot_new(tsh, m, lamda, ttl=ttl)

    fit(tsh, terr, m, A)
    ms = a[1]
    t = a[2]-1

    fit2(tsh, m, t2)

    fit3(tsh, terr, m, t2, A)
    ms3 = a[1]
    t3 = a[2]-1

    fit4(tsh, m, t4)

    fit5(tsh, m, t5)

    t6_0 = t5
    fit6(tsh, terr, m, a)
    ms6 = a[1]

    t6_0 = t
    fit6(tsh, terr, m, a)
    ms7 = a[1]

    t6_0 = t2
    fit6(tsh, terr, m, a)
    ms8 = a[1]

    t6_0 = t4
    fit6(tsh, terr, m, a)
    ms9 = a[1]

    fit10(tsh, m, ms10, t10)

    t6_0 = t5
    fit11(tsh, m, ms11)

    n1 = 5
    n2 = nnn-6

    upg = up[n1:n2]
    vpg = vp[n1:n2]
    arpg = arp[n1:n2]
    umg = um[n1:n2]
    vmg = vm[n1:n2]
    armg = arm[n1:n2]
    tpg = tp[n1:n2]/tb[n1:n2]
    ke = 0.5*avg1d(upg**2+vpg**2)
    pe = 0.5*avg1d((9.8*tpg/n_bar)**2)
    en = ke+pe
    et = en
    tbg = tb[n1:n2]
    kez = 0.5*avg1d(upg**2)
    kev = 0.5*avg1d(vpg**2)
    kevert = 0.5 * avg1d(arpg**2)
    eo = tps_bar/(bo*co)*(9.8/n_bar)**2

    p_rg = p_r[n1:n2]
    rho=p_rg/tbg

    rotary_spec(upg, vpg, acp, ccp, m, lamda, 1, dz*1000)
    up_frac = np.sum(acp)/np.sum(acp+ccp) if lat < 0 else np.sum(ccp)/np.sum(acp+ccp)

    directions_new(upg, vpg, tpg, rho, phi, ut, vt, uq, vq, ui, vi)
    dir = phi*np.rad2deg
    if dir < 0:
        dir=dir+360
    mean_dir=dir

    fac = 9.8/n_bar**2
    uw = fac*w1*ut
    vw = fac*w1*vt
    uw2 = fac*w1*uq
    vw2 = fac*w1*vq
    uwp = -fac*w1p*ut*deltaw
    vwp = -fac*w1p*vt*deltaw
    uw2p = -fac*w1p*uq*deltaw
    vw2p = -fac*w1p*vq*deltaw

    stokes_new(upg, vpg, umg, vmg, phi, n_bar, f_w, df, zm[n1:n2])

    if iw == 0:
        phase_speed_new(umg, vmg, phi, upg, vpg, f, n_bar, f_w, m_bar, k_bar, omega, c_z, c_i, c_x, c_y, df, c_ix, c_iy, c_gx, c_gy, c_xcomp, c_ycomp, u_mean, u_mm, v_mm, theta, dzm)
    if iw == 1:
        phase_speed2(umg, vmg, phi, upg, vpg, f, n_bar, f_w, m_bar, k_bar, omega, c_z, c_i, c_x, c_y, df, c_ix, c_iy, c_gx, c_gy, c_xcomp, c_ycomp, u_mean, u_mm, v_mm, theta)
    return
