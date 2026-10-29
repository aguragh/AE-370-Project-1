# AE370 quadrotor suspended-payload results

This package completes Q1--Q3 using the supplied model and from-scratch RK4. It includes runnable code, raw numerical summaries, figures, and a LaTeX results section. Python requires only NumPy and Matplotlib plus the standard library. No SciPy time integrator is used.

## GitHub updates

- Replace `results.py` with the version here.
- Add `linear_model.py`, `check_linear_model.py`, `results_figures.py`, and `summarize_results.py`.
- Add `data/` and the five `q*.pdf`/`q*.png` figure pairs in `figures/`.
- Keep your current `params.py`, `dynamics.py`, `integrator.py`, `simulate.py`, `verification.py`, and `verification_figures.py`: copies are included unchanged for a self-contained package.
- Merge these instructions into your existing README.
- Put `results_section.tex` and its figures with the report. Replace the old Results outline with `\input{results_section}`. Do not replace your newer Verification section with the older uploaded snippet.

## Commands

```bash
python -m pip install numpy matplotlib
python check_linear_model.py
python results.py
python results_figures.py
python summarize_results.py
```

The sweep uses four worker processes and can take several minutes or longer depending on the computer. The `if __name__ == '__main__'` guard supports Windows. For serial-like resource use, change `max_workers=4` to `max_workers=1` in results.py.

`results.py` reuses existing CSV rows to resume. To recompute after changing the model, gains, thresholds, or simulation implementation, remove the generated data CSV/JSON files first. Otherwise old cached data would be reused. Do not delete the supplied reference outputs unless you intend to regenerate them.

## Output data

- `sweep.csv`: 5,247 runs, 9 masses x 11 lengths x 53 maneuver durations. Angles in radians, lengths in meters, times in seconds, forces in newtons. Boolean columns use 0/1. Failed runs have NaN for full-maneuver metrics; their partial trajectories do not define a completed peak.
- `transitions.json`: each acceptance transition seen on the 0.25 s coarse grid, locally refined to a bracket narrower than 0.005 s at integration step 0.005 s. This does not prove the absence of narrow acceptable intervals between sampled durations.
- `boundary_checks.csv`: both sides of every bracket rerun at steps 0.005 and 0.0025 s. These compare scalar metrics, not nonnested time arrays.
- `corner_checks.csv`: the 12 parameter-extreme cases at both integration steps, reproducing the checks in Section 3.7.
- `window_checks.csv`: four cases with 10 and 20 s observation windows, reproducing the observation-window check.
- `linear_comparison.csv`: full forced hover-linearized model versus the nonlinear model at nominal payload mass and three cable lengths. Relative error uses the nonlinear peak. Saturation is flagged separately.
- `linear_checks.json`: small-perturbation consistency tests of the new forced linearized model.
- `run_log.txt`: progress and measured elapsed time from the execution used to prepare the results.

The original dynamics use unbounded ideal pitch torque, fixed maximum thrust, and no ground model. The results are model/controller-specific and are not hardware-certified safety limits. Cable slack is checked at saved time steps; the model does not simulate a slack cable or retensioning.

## Read the results

Open `results_preview.pdf` for the compiled six-page section. `results_section.tex` is the insertion file for your report; the preview wrapper simply numbers it as Section 4. The figure and table counters in the preview match the uploaded PDF. Your full report will number them automatically.

The main new findings are nominal-payload transition brackets near 7.49, 8.48, and 8.76 s, nonmonotone residual swing for the short cable, and small linear-model peak error in the acceptable-move regime. All 198 refined bracket endpoints retained their classifications at the smaller integration step. These results and their limitations are explained in the section.
