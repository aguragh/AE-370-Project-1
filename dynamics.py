"""Right-hand side f(t, y) of the IVP (eq. 17): reference, controller, mass-matrix solve.

State y = [x, z, theta, phi, xd, zd, thetad, phid]  (eq. 16)
"""
import numpy as np
from params import g


def reference(t, p):
    """Quintic reference x_ref and its first two derivatives (eq. 11)."""
    D, Tm = p["D"], p["Tm"]
    if t >= Tm:
        return D, 0.0, 0.0
    s = t / Tm
    return (D * (10 * s**3 - 15 * s**4 + 6 * s**5),
            D / Tm * (30 * s**2 - 60 * s**3 + 30 * s**4),
            D / Tm**2 * (60 * s - 180 * s**2 + 120 * s**3))


def controller(t, y, p):
    """Cascaded PD, swing-unaware (eqs. 12-15). Returns thrust F and torque tau."""
    x, z, th, _, xd, zd, thd, _ = y
    xr, xrd, xrdd = reference(t, p)
    ax = xrdd + p["kd"] * (xrd - xd) + p["kp"] * (xr - x)                  # (12)
    az = -p["kpz"] * z - p["kdz"] * zd                                     # (13)
    F = min((p["M"] + p["m"]) * np.hypot(ax, g + az), p["Fmax"])           # (14)
    th_cmd = np.arctan2(ax, g + az)                                        # (14)
    tau = p["I"] * (p["kpth"] * (th_cmd - th) - p["kdth"] * thd)          # (15)
    return F, tau


def accelerations(y, F, tau, p):
    """Solve A(phi) [xdd, zdd, phidd]^T = b (eq. 9); thetadd = tau / I (eq. 8)."""
    M, m, l = p["M"], p["m"], p["l"]
    th, ph, phd = y[2], y[3], y[7]
    s, c = np.sin(ph), np.cos(ph)
    A = np.array([[M + m, 0.0, m * l * c],
                  [0.0, M + m, m * l * s],
                  [m * l * c, m * l * s, m * l**2]])
    b = np.array([F * np.sin(th) + m * l * phd**2 * s,
                  F * np.cos(th) - (M + m) * g - m * l * phd**2 * c,
                  -m * g * l * s])
    a1, a2, a3 = np.linalg.solve(A, b)          # no explicit inverse
    return a1, a2, tau / p["I"], a3


def make_rhs(p, ctrl=controller):
    """Return f(t, y). Pass a different ctrl (e.g. constant thrust) for verification tests."""
    def f(t, y):
        F, tau = ctrl(t, y, p)
        xdd, zdd, thdd, phdd = accelerations(y, F, tau, p)
        return np.array([y[4], y[5], y[6], y[7], xdd, zdd, thdd, phdd])
    return f


def make_tension(p, ctrl=controller):
    """Return Tc(t, y) = m[(g + zdd) cos phi - xdd sin phi + l phid^2] (Sec. 1.7)."""
    def Tc(t, y):
        F, tau = ctrl(t, y, p)
        xdd, zdd, _, _ = accelerations(y, F, tau, p)
        return p["m"] * ((g + zdd) * np.cos(y[3]) - xdd * np.sin(y[3]) + p["l"] * y[7]**2)
    return Tc


def thrust(p, ctrl=controller):
    return lambda t, y: ctrl(t, y, p)[0]


def jacobian_hover(p):
    """Analytical hover Jacobian J = df/dy from eqs. (35)-(38)."""
    M, m, l = p["M"], p["m"], p["l"]
    J = np.zeros((8, 8))
    J[0:4, 4:8] = np.eye(4)                                  # kinematic rows
    J[4, 2] = (M + m) * g / M;  J[4, 3] = m * g / M          # (35)
    J[5, 1] = -p["kpz"];        J[5, 5] = -p["kdz"]          # (36)
    J[6, 0] = -p["kpth"] * p["kp"] / g                       # (37)
    J[6, 4] = -p["kpth"] * p["kd"] / g
    J[6, 2] = -p["kpth"];       J[6, 6] = -p["kdth"]
    J[7, 2] = J[7, 3] = -(M + m) * g / (M * l)               # (38)
    return J
