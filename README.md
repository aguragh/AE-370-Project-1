# Quadrotor with a cable-suspended payload (AE 370 Project 1)

Python 3 with `numpy`. RK4 is implemented from scratch in `integrator.py`; no outside time integrators are used.

| File | Contents |
|---|---|
| `params.py` | Parameter values (Table 1) and controller gains (eq. 21) |
| `dynamics.py` | Right-hand side f(t, y): reference (eq. 11), controller (eqs. 12-15), mass-matrix solve (eq. 9), cable tension, hover Jacobian (eqs. 35-38) |
| `integrator.py` | RK4 step (eqs. 23-27), segmented time grid (Sec. 2.6), time marching with tension check (Algorithm 1), stability function and step limit (eq. 34, Sec. 2.5) |
| `simulate.py` | One maneuver plus quantities of interest and acceptance criteria (Sec. 1.8) |
| `verification.py` | Jacobian check, eigenvalues, stability limit, self-convergence, exact tests (Sec. 3) |
| `results.py` | Parameter sweeps for Q1-Q3 |

```
python verification.py
python results.py
```
