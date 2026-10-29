"""Explicit RK4 time integrator (Sec. 2.2, Algorithm 1). Written from scratch."""
import numpy as np


def rk4_step(f, t, y, h):
    """One RK4 step from t_k to t_k+1 (eqs. 23-27)."""
    k1 = f(t, y)
    k2 = f(t + h / 2, y + h / 2 * k1)
    k3 = f(t + h / 2, y + h / 2 * k2)
    k4 = f(t + h, y + h * k3)
    return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def time_grid(segments, dt):
    """Uniform steps within each segment, joins exactly on grid points (Sec. 2.6).
    segments = [Tm, Tobs] -> N_i = ceil(T_i / dt), h_i = T_i / N_i."""
    t, t0 = [0.0], 0.0
    for T in segments:
        N = int(np.ceil(T / dt - 1e-12))
        t.extend(t0 + T * np.arange(1, N + 1) / N)
        t0 += T
    return np.array(t)


def integrate(f, y0, t, check=None):
    """March y' = f(t, y) over the grid t. Returns (t, Y, ok).

    check(t, y) -> bool is evaluated after every step (Algorithm 1, lines 8-9);
    the run stops and is flagged (ok = False) if it fails or the state is nonfinite."""
    Y = np.empty((len(t), len(y0)))
    Y[0] = y0
    for k in range(len(t) - 1):
        Y[k + 1] = rk4_step(f, t[k], Y[k], t[k + 1] - t[k])
        if not np.all(np.isfinite(Y[k + 1])) or (check and not check(t[k + 1], Y[k + 1])):
            return t[:k + 2], Y[:k + 2], False
    return t, Y, True


def R(w):
    """RK4 stability function (eq. 34)."""
    return 1 + w + w**2 / 2 + w**3 / 6 + w**4 / 24


def stable_step_limit(eigs, hmax=10.0, tol=1e-8):
    """Largest h with |R(h lambda_j)| < 1 for all eigenvalues (Sec. 2.5), by bisection."""
    stable = lambda h: np.all(np.abs(R(h * np.asarray(eigs))) < 1)
    lo, hi = 0.0, hmax
    while hi - lo > tol:
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if stable(mid) else (lo, mid)
    return lo
