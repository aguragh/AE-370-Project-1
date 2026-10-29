"""Reproduce verification figures using the repository's unmodified solver.
Run: python verification_figures.py (requires NumPy and Matplotlib).
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from params import P, g
from dynamics import make_rhs, make_tension
from integrator import time_grid, integrate
from simulate import run

OUT = Path(__file__).resolve().parent / 'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .25})
p = dict(P)
M, m, l = p['M'], p['m'], p['l']

def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / (name + '.pdf'), bbox_inches='tight')
    fig.savefig(OUT / (name + '.png'), dpi=180, bbox_inches='tight')
    plt.close(fig)

def constant(F):
    return lambda t, y, p: (F, 0.)

def open_loop(y0, F, T, h):
    ctrl = constant(F)
    t, Y, ok = integrate(make_rhs(p, ctrl), y0, time_grid([T], h))
    tension = np.array([make_tension(p, ctrl)(ti, yi) for ti, yi in zip(t, Y)])
    assert ok and tension.min() > 0
    return t, Y, float(tension.min())

def freefall_exact(t, rate=2.):
    """COM ballistic trajectory and constant relative rotation, y0[7]=rate."""
    a = m * l / (M + m)
    ph = rate * t
    Y = np.zeros((len(t), 8))
    Y[:, 0] = a * rate * t - a * np.sin(ph)
    Y[:, 1] = -a - .5 * g * t**2 + a * np.cos(ph)
    Y[:, 3] = ph
    Y[:, 4] = a * rate * (1 - np.cos(ph))
    Y[:, 5] = -g * t - a * rate * np.sin(ph)
    Y[:, 7] = rate
    return Y

def main():
    data = {}
    # Baseline hover is a separate equilibrium check with fixed reference.
    ph = dict(p, D=0.)
    t, Y, ok = integrate(make_rhs(ph), np.zeros(8), time_grid([20.], .01))
    assert ok
    data['hover_max_abs_state'] = float(np.abs(Y).max())
    # Nonlinear conservation check at a large initial angle.
    y0 = np.zeros(8); y0[3] = np.pi / 3
    tm, Ym, tc_m = open_loop(y0, (M + m) * g, 20., .001)
    px = (M + m) * Ym[:, 4] + m * l * Ym[:, 7] * np.cos(Ym[:, 3])
    pm = px - px[0]
    data['momentum_max_drift_kg_m_s'] = float(abs(pm).max())
    data['momentum_min_tension_N'] = tc_m
    # Exact free fall tests all translating state components, not only constant rate.
    y0 = np.zeros(8); y0[7] = 2.
    tf, Yf, tc_f = open_loop(y0, 0., 3., .001)
    exact = freefall_exact(tf)
    x, z, th, phi, vx, vz, om, rate = Yf.T
    energy = (.5*M*(vx**2 + vz**2) + .5*m*((vx+l*rate*np.cos(phi))**2
              + (vz+l*rate*np.sin(phi))**2) + .5*p['I']*om**2
              + (M+m)*g*z-m*g*l*np.cos(phi))
    de = energy - energy[0]
    data['freefall_max_energy_drift_J'] = float(abs(de).max())
    data['freefall_max_rate_error_rad_s'] = float(abs(rate - 2.).max())
    data['freefall_max_position_error_m'] = float(np.linalg.norm(Yf[:,:2]-exact[:,:2],axis=1).max())
    data['freefall_min_tension_N'] = tc_f
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
    ax[0].plot(tm, pm); ax[0].set(xlabel='Time (s)', ylabel='Momentum drift (kg m/s)', title='(a) Nonlinear momentum test')
    ax[1].plot(tf, de); ax[1].set(xlabel='Time (s)', ylabel='Energy drift (J)', title='(b) Free-fall energy test')
    ax[2].plot(tf, Yf[:,0], label='RK4 x'); ax[2].plot(tf[::120], exact[::120,0], 'o', ms=3, label='Exact x')
    ax[2].plot(tf, Yf[:,1], label='RK4 z'); ax[2].plot(tf[::120], exact[::120,1], 'x', ms=4, label='Exact z')
    ax[2].set(xlabel='Time (s)', ylabel='Vehicle position (m)', title='(c) Exact free-fall trajectory'); ax[2].legend(fontsize=8)
    save(fig, 'verification_invariants')
    # Approximation check: the linear period is exact only as amplitude -> 0.
    omega = np.sqrt(g*(M+m)/(M*l)); Tlin = 2*np.pi/omega
    periods = []; small = None
    for deg in [1, 10, 30, 60]:
        y0 = np.zeros(8); y0[3] = np.deg2rad(deg)
        t, Y, tc = open_loop(y0, (M+m)*g, 30., .001)
        idx = np.where((Y[:-1,3]<0) & (Y[1:,3]>=0))[0]
        crossings = t[idx] - Y[idx,3]*(t[idx+1]-t[idx])/(Y[idx+1,3]-Y[idx,3])
        period = float(np.diff(crossings).mean())
        periods.append([deg, period, 100*(period/Tlin-1), tc])
        if deg == 1: small = t, Y[:,3]
    data['small_angle_period_s'] = float(Tlin)
    data['period_rows_angle_deg_period_s_error_percent_min_tension_N'] = periods
    fig, ax = plt.subplots(1,2,figsize=(9,3.4))
    t, phi = small; mask = t<=8
    ax[0].plot(t[mask], np.rad2deg(phi[mask]), label='Nonlinear RK4')
    ax[0].plot(t[mask], np.cos(omega*t[mask]), '--', label='Small-angle approximation')
    ax[0].set(xlabel='Time (s)', ylabel='Swing angle (deg)', title='(a) Initial swing = 1 degree'); ax[0].legend(fontsize=8)
    ax[1].plot(np.array(periods)[:,0], np.array(periods)[:,2], 'o-')
    ax[1].set(xlabel='Initial swing angle (deg)', ylabel='Period difference (%)', title='(b) Departure from small-angle period')
    save(fig,'verification_small_angle')
    # Common-time self convergence; nesting is verified for this baseline only.
    hs = [.04,.02,.01,.005]; runs = [run(p,h) for h in hs]
    errs = []
    for a,b in zip(runs[:-1],runs[1:]):
        assert a['completed'] and b['completed']
        assert np.allclose(a['t'], b['t'][::2], atol=1e-12, rtol=0)
        errs.append(float(abs(a['Y'][:,3]-b['Y'][::2,3]).max()))
    data['self_steps_s'] = hs[:-1]; data['self_errors_rad'] = errs
    data['self_orders'] = np.log2(np.array(errs[:-1])/errs[1:]).tolist()
    data['peak_difference_deg_001_0005'] = float(np.rad2deg(abs(runs[2]['phi_peak']-runs[3]['phi_peak'])))
    # True error in vehicle position for exact nonlinear free fall.
    he = [.06,.03,.015,.0075]; ee = []
    y0 = np.zeros(8); y0[7] = 2.
    for h in he:
        t,Y,_ = open_loop(y0,0.,3.,h)
        ee.append(float(np.linalg.norm(Y[:,:2]-freefall_exact(t)[:,:2],axis=1).max()))
    data['exact_steps_s'] = he; data['exact_position_errors_m'] = ee
    data['exact_orders'] = np.log2(np.array(ee[:-1])/ee[1:]).tolist()
    fig, ax = plt.subplots(1,2,figsize=(9,3.4))
    for a, h, e, title, label in [(ax[0],hs[:-1],errs,'(a) Controlled maneuver','Common-time swing difference (rad)'),
                                 (ax[1],he,ee,'(b) Exact free-fall solution','Maximum position error (m)')]:
        h=np.array(h); e=np.array(e)
        a.loglog(h,e,'o-',label='Measured')
        a.loglog(h,e[0]*(h/h[0])**4,'--',label='Fourth-order guide')
        a.set(xlabel='Time step (s)',ylabel=label,title=title); a.legend(fontsize=8)
        a.set_xticks(h, [f'{value:g}' for value in h], fontsize=9)
        a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    save(fig,'verification_convergence')
    (OUT.parent/'verification_metrics.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))

if __name__ == '__main__':
    main()
