"""Check the forcing terms in the full hover linearization independently.
Run python check_linear_model.py. Saves data/linear_checks.json.
"""
from pathlib import Path
import json
import numpy as np
from params import P
from dynamics import make_rhs
from linear_model import make_linear_rhs,run_linear
from simulate import run

v=np.array([.2,-.1,.03,-.04,.1,-.05,.02,.03])
checks=[]
for e in [1e-2,1e-3,1e-4]:
    p=dict(P,D=e*20)
    err=float(np.max(abs(make_rhs(p)(3,e*v)-make_linear_rhs(p)(3,e*v))))
    checks.append(dict(perturbation=e,rhs_error=err))
# First-order expansion must leave second-order or smaller residuals.
assert checks[0]['rhs_error']/checks[1]['rhs_error']>90
assert checks[1]['rhs_error']/checks[2]['rhs_error']>90
trajectory=[]
for D in [.2,.02]:
    p=dict(P,D=D);r=run(p);t,Y=run_linear(p)
    assert r['completed'] and np.allclose(t,r['t'])
    trajectory.append(dict(D=D,swing_error_rad=float(np.max(abs(r['Y'][:,3]-Y[:,3])))))
assert trajectory[0]['swing_error_rad']<2e-9
result=dict(rhs_checks=checks,trajectory_checks=trajectory)
path=Path(__file__).resolve().parent/'data'/'linear_checks.json'
path.write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
