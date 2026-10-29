"""Parameter values (Table 1) and controller gains (eq. 21)."""
import numpy as np

g = 9.81

# Baseline case (Sec. 2.1): m = 0.5 kg, l = 1.5 m, D = 20 m, Tm = 10 s
P = dict(M=1.5, I=0.02, m=0.5, l=1.5, D=20.0, Tm=10.0, Tobs=10.0,
         Fmax=2 * (1.5 + 0.5) * g,   # fixed at 39.24 N using nominal m0 = 0.5 kg
         kp=0.5, kd=1.4,             # horizontal position
         kpz=4.0, kdz=4.0,           # altitude
         kpth=100.0, kdth=20.0)      # attitude

# Sweep ranges (Table 1)
M_PAYLOAD = np.linspace(0.2, 1.0, 9)    # kg
L_CABLE = np.linspace(0.5, 3.0, 11)     # m
T_MOVE = np.arange(2.0, 15.01, 0.25)    # s (refine near transitions later)
PHI_MAX = np.deg2rad(15.0)
