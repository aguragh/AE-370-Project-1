"""Full first-order closed-loop hover linearization, with the same reference."""
import numpy as np
from params import g
from dynamics import jacobian_hover, reference
from integrator import time_grid, integrate


def make_linear_rhs(p):
    J = jacobian_hover(p)
    def f(t, y):
        xr, vr, ar = reference(t, p)
        dy = J @ y
        dy[6] += p['kpth'] / g * (ar + p['kd'] * vr + p['kp'] * xr)
        return dy
    return f


def run_linear(p, dt=.01):
    t, Y, ok = integrate(make_linear_rhs(p), np.zeros(8), time_grid([p['Tm'],p['Tobs']],dt))
    assert ok
    return t, Y
