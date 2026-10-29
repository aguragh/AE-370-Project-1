# AE370: Quadrotor with a Cable-Suspended Payload

This project studies how maneuver duration, cable length, and payload
mass affect the swing of a payload carried by a planar quadrotor.

The nonlinear equations are integrated using a fourth-order
Runge–Kutta method implemented from scratch. The study examines:

1. The shortest resolved maneuver durations satisfying the swing,
   cable-tension, altitude, and final tracking requirements.
2. Residual payload swing after the commanded move ends.
3. Agreement between the nonlinear model and its full closed-loop
   hover linearization.

## Requirements

The code uses Python, NumPy, and Matplotlib.

Install the required packages from the repository directory:

```bash
python -m pip install -r requirements.txt
```

Alternatively:

```bash
python -m pip install numpy matplotlib
```

## Repository contents

| File or folder | Purpose |
|---|---|
| `params.py` | Physical parameters, controller gains, and sweep ranges |
| `dynamics.py` | Nonlinear equations, reference trajectory, controller, and cable tension |
| `integrator.py` | From-scratch RK4 method and time-grid construction |
| `simulate.py` | Runs one nonlinear maneuver and evaluates its performance |
| `verification.py` | Jacobian, stability, conservation, and convergence checks |
| `verification_figures.py` | Generates verification figures and numerical summaries |
| `linear_model.py` | Full forced closed-loop hover linearization |
| `check_linear_model.py` | Checks consistency of the linearized and nonlinear models |
| `results.py` | Parameter sweep, boundary refinement, and sensitivity checks |
| `results_figures.py` | Generates results figures from saved data |
| `summarize_results.py` | Prints numerical summaries used in the report |
| `requirements.txt` | Required Python packages |
| `data/` | Saved sweep results, boundary checks, and linear-model comparisons |
| `figures/` | Generated verification and results figures |

## Run the verification

From the repository directory, run:

```bash
python verification.py
python verification_figures.py
python check_linear_model.py
```

These scripts check hover equilibrium, the hover Jacobian, momentum
and energy conservation, the small-angle swing limit, exact free-fall
motion, RK4 convergence, and the forced linearized model.

Verification figures are written to `figures/`.
Numerical summaries are written to `verification_metrics.json`
and `data/linear_checks.json`.

Small differences in floating-point residuals between computers
are expected.

## View the results without rerunning simulations

The generated data are included in `data/`.

To recreate the five results figures and print the numerical summaries:

```bash
python results_figures.py
python summarize_results.py
```

The figures are saved in `figures/` as PDF and PNG files.

## Run the full results study

```bash
python results.py
python results_figures.py
python summarize_results.py
```

The coarse sweep contains 5,247 nonlinear simulations:

- Nine payload masses from 0.2 to 1.0 kg.
- Eleven cable lengths from 0.5 to 3.0 m.
- Fifty-three maneuver durations from 2 to 15 s.

The script also refines detected acceptance transitions, checks their
sensitivity to the integration step, reproduces the parameter-extreme
and observation-window checks, and compares nonlinear predictions
with the linearized model.

Four worker processes are used. Execution time depends on the computer.
The original complete data-generation run took approximately
7.6 minutes in the execution environment.

### Cached data and fresh runs

`results.py` reuses completed rows from existing output CSV files.
This allows an interrupted study to resume.

For a completely fresh run, or after changing the model, parameters,
controller, or acceptance criteria, remove these generated files
from `data/` before running `results.py`:

- `sweep.csv`
- `transitions.json`
- `boundary_checks.csv`
- `corner_checks.csv`
- `window_checks.csv`
- `linear_comparison.csv`

Keep a copy of the original outputs if you want to compare results.
Do not mix cached outputs from different model versions.

## Main results files

| Output | Contents |
|---|---|
| `sweep.csv` | Coarse parameter sweep and acceptance classifications |
| `transitions.json` | Refined fail/pass brackets in maneuver duration |
| `boundary_checks.csv` | Comparisons at integration steps of 0.005 and 0.0025 s |
| `corner_checks.csv` | Parameter-extreme comparisons at two integration steps |
| `window_checks.csv` | Comparisons using 10 and 20 s observation intervals |
| `linear_comparison.csv` | Linear and nonlinear peak-swing predictions |
| `linear_checks.json` | Small-perturbation checks of the linearized model |

In the nonlinear sweep CSV files, angles are in radians, times in
seconds, lengths in meters, and forces in newtons. Boolean fields
use 0 and 1. The linear-comparison file labels angle columns in degrees.

Runs stopped by the cable-tension check have `NaN` for full-maneuver
metrics because the simulation did not complete.

## Interpretation and limitations

Acceptable maneuvers satisfy all specified conditions, not only the
15-degree swing limit.

Reported maneuver-time boundaries are refined transitions detected
from a finite parameter grid. They are not globally optimal maneuver
times or hardware-certified safety limits.

The model assumes planar motion, a massless inextensible cable,
ideal pitch torque, and no aerodynamic drag or ground contact.
Cable tension is checked at saved time points; motion after
detected cable slack is not simulated.
