"""Parameter sweep for Q1-Q2 (Sec. 1.8). Start coarse, then refine near transitions."""
import numpy as np
from params import P, T_MOVE, L_CABLE
from simulate import run

p = dict(P)
results = {}
for l in L_CABLE[::5]:                     # coarse first: l = 0.5, 1.75, 3.0 m
    rows = []
    for Tm in T_MOVE:
        r = run(dict(p, l=l, Tm=Tm))
        rows.append((Tm, np.rad2deg(r["phi_peak"]), np.rad2deg(r["phi_res"]), r["acceptable"]))
    rows = np.array(rows)
    results[l] = rows
    ok = rows[rows[:, 3] == 1, 0]
    print(f"l = {l:.2f} m: Tm* = {ok.min() if ok.size else 'none in [2, 15]'} s")

np.savez("sweep_m0.5.npz", **{f"l_{l:.2f}": v for l, v in results.items()})
