"""Run one maneuver and compute the quantities of interest (Sec. 1.8)."""
import numpy as np
from dynamics import make_rhs, make_tension, thrust
from integrator import time_grid, integrate


def run(p, dt=0.01):
    f, Tc, Fn = make_rhs(p), make_tension(p), thrust(p)
    t = time_grid([p["Tm"], p["Tobs"]], dt)
    t, Y, ok = integrate(f, np.zeros(8), t, check=lambda tk, yk: Tc(tk, yk) > 0)
    tens = np.array([Tc(tk, yk) for tk, yk in zip(t, Y)])
    F = np.array([Fn(tk, yk) for tk, yk in zip(t, Y)])
    phi = np.abs(Y[:, 3])
    last = t >= t[-1] - 1.0                                  # final 1 s
    pos_err = np.hypot(Y[last, 0] - p["D"], Y[last, 1]).max()
    speed = np.hypot(Y[last, 4], Y[last, 5]).max()
    out = dict(t=t, Y=Y, completed=ok,
               phi_peak=phi.max(),                           # eq. (19)
               phi_res=phi[t >= p["Tm"]].max() if ok else np.nan,
               Tc_min=tens.min(), F_max=F.max(), saturated=F.max() >= p["Fmax"] - 1e-9,
               z_max=np.abs(Y[:, 1]).max(), pos_err=pos_err, speed=speed)
    out["acceptable"] = (ok and out["phi_peak"] <= np.deg2rad(15) and out["z_max"] <= 0.5
                         and pos_err <= 0.2 and speed <= 0.2)
    return out
