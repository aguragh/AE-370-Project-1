"""Section 2.5-2.6 and Section 3 checks: Jacobian, stability limit, exact tests, convergence."""
import numpy as np
from params import P, g, M_PAYLOAD, L_CABLE
from dynamics import make_rhs, jacobian_hover
from integrator import time_grid, integrate, stable_step_limit
from simulate import run

p = dict(P)
M, m, l = p["M"], p["m"], p["l"]

# ---- Sec. 2.5: analytical Jacobian vs centered differences, eigenvalues, stability limit ----
p_h = dict(p, D=0.0)                       # fixed hover reference at the origin (t > Tm)
f = make_rhs(p_h)
J = jacobian_hover(p_h)
eps, Jfd = 1e-6, np.zeros((8, 8))
for j in range(8):
    e = np.zeros(8); e[j] = eps
    Jfd[:, j] = (f(20.0, e) - f(20.0, -e)) / (2 * eps)
eigs = np.linalg.eigvals(J)
print(f"max |J - J_fd|          = {np.abs(J - Jfd).max():.1e}")
print("eigenvalues             =", np.round(np.sort_complex(eigs), 4))
print(f"RK4 stable step limit   = {stable_step_limit(eigs):.4f} s")

worst_re, worst_h = -np.inf, np.inf
for mm in M_PAYLOAD:
    for ll in L_CABLE:
        ev = np.linalg.eigvals(jacobian_hover(dict(p, m=mm, l=ll)))
        worst_re = max(worst_re, ev.real.max())
        worst_h = min(worst_h, stable_step_limit(ev))
print(f"grid: max Re(lambda)    = {worst_re:.4f},  min step limit = {worst_h:.4f} s")

# ---- Sec. 2.6: self-convergence of phi on common grid points ----
dts = [0.04, 0.02, 0.01, 0.005]
runs = [run(p, dt) for dt in dts]
E = [np.abs(a["Y"][:, 3] - b["Y"][::2, 3]).max() for a, b in zip(runs[:-1], runs[1:])]
for dt, e in zip(dts, E):
    print(f"dt = {dt:<6} max|phi_h - phi_h/2| = {e:.4e} rad")
print("observed order p        =", np.round(np.log2(np.array(E[:-1]) / E[1:]), 3))
b = runs[2]
print(f"baseline dt=0.01: peak = {np.rad2deg(b['phi_peak']):.3f} deg, "
      f"Tc_min = {b['Tc_min']:.3f} N, F_max = {b['F_max']:.3f} N, "
      f"pos_err = {b['pos_err']:.5f} m, speed = {b['speed']:.5f} m/s")
print(f"peak change 0.01->0.005 = {np.rad2deg(abs(runs[3]['phi_peak'] - b['phi_peak'])):.5f} deg")

# ---- Sec. 3: exact special cases (open-loop, constant thrust, zero torque) ----
const = lambda F: (lambda t, y, p: (F, 0.0))
grid = lambda T, dt=1e-3: time_grid([T], dt)

y0 = np.zeros(8); y0[3] = np.deg2rad(60)                                   # Test 1
_, Y, _ = integrate(make_rhs(p, const((M + m) * g)), y0, grid(20.0))
px = (M + m) * Y[:, 4] + m * l * Y[:, 7] * np.cos(Y[:, 3])
print(f"Test 1 momentum: max |px - px0| = {np.abs(px - px[0]).max():.2e}")

T_exact = 2 * np.pi / np.sqrt(g * (M + m) / (M * l))                        # Test 2
for deg in [1, 10, 30, 60]:
    y0 = np.zeros(8); y0[3] = np.deg2rad(deg)
    t, Y, _ = integrate(make_rhs(p, const((M + m) * g)), y0, grid(30.0))
    up = np.where((Y[:-1, 3] < 0) & (Y[1:, 3] >= 0))[0]
    tc = t[up] - Y[up, 3] * (t[up + 1] - t[up]) / (Y[up + 1, 3] - Y[up, 3])   # interpolated crossings
    T_sim = np.mean(np.diff(tc))
    print(f"Test 2 phi0 = {deg:>2} deg: period {T_sim:.5f} s vs linear {T_exact:.5f} s "
          f"({100 * (T_sim / T_exact - 1):+.3f}%)")

y0 = np.zeros(8); y0[7] = 2.0                                              # Test 3
_, Y, _ = integrate(make_rhs(p, const(0.0)), y0, grid(3.0))
ph, phd, xd, zd, z = Y[:, 3], Y[:, 7], Y[:, 4], Y[:, 5], Y[:, 1]
E_tot = (0.5 * M * (xd**2 + zd**2)
         + 0.5 * m * ((xd + l * phd * np.cos(ph))**2 + (zd + l * phd * np.sin(ph))**2)
         + M * g * z + m * g * (z - l * np.cos(ph)))
print(f"Test 3 free fall: max |phid - phid0| = {np.abs(phd - phd[0]).max():.1e}, "
      f"max |E - E0| = {np.abs(E_tot - E_tot[0]).max():.1e}")
