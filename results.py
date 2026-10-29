"""Reproduce Q1--Q3 data. NumPy only; RK4 and dynamics are unchanged.
Run python results.py. Resume by rerunning; completed CSV rows are reused.
Run python results_figures.py afterward to generate plots.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import time
import numpy as np
from params import P, M_PAYLOAD, L_CABLE, T_MOVE
from simulate import run
from linear_model import run_linear

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
DATA.mkdir(exist_ok=True)
FIELDS = ['m','l','Tm','dt','Tobs','completed','acceptable','saturated',
          'phi_peak','phi_res','Tc_min','F_max','z_max','pos_err','speed','runtime_s']


def case(args):
    m,l,Tm,dt,obs = args
    p = dict(P,m=float(m),l=float(l),Tm=float(Tm),Tobs=float(obs))
    start=time.perf_counter(); r=run(p,dt)
    row=dict(m=m,l=l,Tm=Tm,dt=dt,Tobs=obs,runtime_s=time.perf_counter()-start)
    for k in FIELDS[5:-1]: row[k]=int(r[k]) if isinstance(r[k],(bool,np.bool_)) else float(r[k])
    # A terminated trajectory does not define full-maneuver peaks or arrival metrics.
    if not r['completed']:
        for k in ['phi_peak','phi_res','z_max','pos_err','speed']: row[k]=float('nan')
    return row


def read(name):
    path=DATA/name
    if not path.exists(): return []
    with path.open() as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]


def batch(name, jobs):
    old=read(name)
    key=lambda r:tuple(round(r[k],10) for k in FIELDS[:5])
    done={key(r) for r in old}
    jobs=[j for j in jobs if tuple(round(float(x),10) for x in j) not in done]
    path=DATA/name
    with path.open('a',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS)
        if not old:w.writeheader()
        with ProcessPoolExecutor(max_workers=4) as pool:
            for i,r in enumerate(pool.map(case,jobs,chunksize=4)):
                w.writerow(r); f.flush();old.append(r)
                if (i+1)%200==0: print(name,i+1,'/',len(jobs),flush=True)
    return old


def refine_group(args):
    m,l,rows=args
    rows=sorted(rows,key=lambda r:r['Tm'])
    transitions=[]
    for a,b in zip(rows[:-1],rows[1:]):
        if bool(a['acceptable'])==bool(b['acceptable']):continue
        lo,hi=a['Tm'],b['Tm']; left=bool(a['acceptable'])
        # Local transition refinement only: no global monotonicity assumption.
        while hi-lo>.005:
            mid=(lo+hi)/2
            r=case((m,l,mid,.005,10.))
            if bool(r['acceptable'])==left:lo=mid
            else:hi=mid
        ra=case((m,l,lo,.005,10.));rb=case((m,l,hi,.005,10.))
        transitions.append(dict(m=m,l=l,lo=lo,hi=hi,left_acceptable=int(ra['acceptable']),
                                right_acceptable=int(rb['acceptable']),
                                phi_lo_deg=float(np.rad2deg(ra['phi_peak'])),
                                phi_hi_deg=float(np.rad2deg(rb['phi_peak']))))
    return transitions


def main():
    start=time.perf_counter()
    # 9 payload masses x 11 cable lengths x 53 maneuver durations.
    rows=batch('sweep.csv',[(m,l,T,.01,10.) for m in M_PAYLOAD for l in L_CABLE for T in T_MOVE])
    print('Sweep complete:',len(rows),flush=True)
    groups=[(m,l,[r for r in rows if abs(r['m']-m)<1e-8 and abs(r['l']-l)<1e-8])
            for m in M_PAYLOAD for l in L_CABLE]
    transitions=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for i,t in enumerate(pool.map(refine_group,groups)):
            transitions.extend(t)
            if (i+1)%10==0:print('Boundary groups:',i+1,'/99',flush=True)
    (DATA/'transitions.json').write_text(json.dumps(transitions,indent=2))
    # Validate both sides of every refined transition at a further halved step.
    jobs=[]
    for r in transitions:
        for T in [r['lo'],r['hi']]:
            for h in [.005,.0025]:jobs.append((r['m'],r['l'],T,h,10.))
    batch('boundary_checks.csv',jobs)
    # Reproduce Section 3.7 corner checks and observation-window comparisons.
    batch('corner_checks.csv',[(m,l,T,h,10.) for m in [.2,1.] for l in [.5,3.]
                              for T in [2.,6.,15.] for h in [.01,.005]])
    batch('window_checks.csv',[(.5,l,T,.01,obs) for l,T in [(.5,7.5),(1.75,8.5),(3.,9.),(1.5,10.)]
                              for obs in [10.,20.]])
    # Full linearized-model comparison at nominal payload, three cable lengths.
    out=[]
    for l in [.5,1.75,3.]:
        for T in T_MOVE:
            _,Y=run_linear(dict(P,l=l,Tm=T))
            r=next(r for r in rows if abs(r['m']-.5)<1e-8 and abs(r['l']-l)<1e-8 and r['Tm']==T)
            lin=float(np.max(np.abs(Y[:,3])))
            valid=bool(r['completed']); err=abs(lin-r['phi_peak']) if valid else float('nan')
            out.append(dict(l=l,Tm=T,nonlinear_deg=float(np.rad2deg(r['phi_peak'])),linear_deg=float(np.rad2deg(lin)),
                            abs_error_deg=float(np.rad2deg(err)),relative_percent=100*err/r['phi_peak'] if valid and r['phi_peak']>1e-6 else float('nan'),
                            completed=int(valid),saturated=int(r['saturated'])))
    with (DATA/'linear_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
    print('Data complete. This invocation wall time:',time.perf_counter()-start,flush=True)

if __name__=='__main__':main()
