"""Generate report figures from saved CSV/JSON data; no simulations rerun."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parent
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2})

def read(name):
    with (ROOT/'data'/name).open() as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]

def save(fig,name):
    fig.tight_layout()
    for ext in ['pdf','png']:fig.savefig(FIG/(name+'.'+ext),dpi=200,bbox_inches='tight')
    plt.close(fig)

def group(rows,m,l):return sorted([r for r in rows if abs(r['m']-m)<1e-8 and abs(r['l']-l)<1e-8],key=lambda r:r['Tm'])

def main():
    rows=read('sweep.csv')
    trans=json.loads((ROOT/'data'/'transitions.json').read_text())
    colors=['tab:blue','tab:orange','tab:green']; ls=[.5,1.75,3.]
    fig,ax=plt.subplots(1,2,figsize=(10,3.9))
    for l,c in zip(ls,colors):
        rr=group(rows,.5,l);t=np.array([r['Tm'] for r in rr]);p=np.rad2deg([r['phi_peak'] for r in rr])
        for a in ax:a.plot(t,p,color=c,label=f'l = {l:g} m')
        ok=np.array([bool(r['acceptable']) for r in rr]);ax[1].plot(t[ok],p[ok],'.',color=c,ms=4)
        bad=np.array([not bool(r['completed']) for r in rr]);ax[0].plot(t[bad],np.zeros(sum(bad)),'x',color=c,ms=5)
    for a in ax:
        a.axhline(15,color='k',ls='--',lw=1,label='15-degree limit');a.set(xlabel='Maneuver duration (s)',ylabel='Peak swing (deg)')
    ax[0].set(title='(a) Full range; crosses mark stopped runs',xlim=(2,15));ax[0].legend(fontsize=8)
    ax[1].set(title='(b) Near the swing limit; dots are acceptable',xlim=(6,11),ylim=(5,35))
    save(fig,'q1_peak_swing')
    # Map nominal mass: failed runs are masked, not treated as full-trajectory peaks.
    rr=[r for r in rows if abs(r['m']-.5)<1e-8]
    cables=sorted(set(r['l'] for r in rr));times=sorted(set(r['Tm'] for r in rr))
    peak=np.array([[np.rad2deg(r['phi_peak']) for r in group(rows,.5,l)] for l in cables])
    status=np.array([[0 if not r['completed'] else (2 if r['acceptable'] else 1) for r in group(rows,.5,l)] for l in cables])
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    mesh=ax[0].pcolormesh(times,cables,np.ma.masked_invalid(peak),shading='nearest',rasterized=True,vmin=0,vmax=60,cmap='viridis')
    ax[0].set_facecolor('.8');fig.colorbar(mesh,ax=ax[0],label='Peak swing (deg)',extend='max')
    ax[0].contour(times,cables,np.ma.masked_invalid(peak),levels=[15],colors='white',linewidths=1.5)
    ax[0].set(title='(a) Peak swing; white contour = 15 degrees')
    ax[1].pcolormesh(times,cables,status,shading='nearest',rasterized=True,cmap=ListedColormap(['.7','#e6a07a','#74b69b']),vmin=0,vmax=2)
    ax[1].set(title='(b) All acceptance conditions')
    ax[1].legend(handles=[Patch(color=c,label=s) for c,s in [('.7','Stopped'),('#e6a07a','Unacceptable'),('#74b69b','Acceptable')]],fontsize=8,loc='upper right')
    for a in ax:a.set(xlabel='Maneuver duration (s)',ylabel='Cable length (m)',xlim=(2,15),ylim=(.5,3))
    save(fig,'q1_acceptance_map')
    fig,ax=plt.subplots(figsize=(6.5,4))
    for l,c in zip(ls,colors):
        rr=sorted([r for r in trans if abs(r['l']-l)<1e-8 and r['right_acceptable'] and not r['left_acceptable']],key=lambda r:r['m'])
        ax.errorbar([r['m']/1.5 for r in rr],[(r['lo']+r['hi'])/2 for r in rr],yerr=[(r['hi']-r['lo'])/2 for r in rr],fmt='o-',color=c,label=f'l = {l:g} m',capsize=3)
    ax.set(xlabel='Payload-to-vehicle mass ratio m/M',ylabel='First resolved acceptable duration (s)');ax.legend()
    save(fig,'q1_mass_effect')
    fig,ax=plt.subplots(1,2,figsize=(10,3.8))
    for l,c in zip(ls,colors):
        rr=group(rows,.5,l);t=np.array([r['Tm'] for r in rr]);res=np.rad2deg([r['phi_res'] for r in rr])
        ax[0].plot(t,res,color=c,label=f'l = {l:g} m');ax[1].plot(t,res,color=c)
        ok=np.array([bool(r['acceptable']) for r in rr]);ax[1].plot(t[ok],res[ok],'.',color=c,ms=4)
    for a in ax:a.set(xlabel='Maneuver duration (s)',ylabel='Residual peak swing (deg)')
    ax[0].set(title='(a) All completed runs');ax[0].legend(fontsize=8)
    ax[1].set(title='(b) Longer moves; dots are acceptable',xlim=(7,15),ylim=(0,10))
    save(fig,'q2_residual_swing')
    lin=read('linear_comparison.csv')
    fig,ax=plt.subplots(1,3,figsize=(12,3.6))
    for l,c in zip(ls,colors):
        rr=[r for r in lin if r['l']==l];t=np.array([r['Tm'] for r in rr])
        ax[0].plot(t,[r['nonlinear_deg'] for r in rr],color=c,label=f'l = {l:g} m')
        ax[0].plot(t,[r['linear_deg'] if r['completed'] else np.nan for r in rr],color=c,ls='--')
        for a,k in [(ax[1],'abs_error_deg'),(ax[2],'relative_percent')]:
            a.plot(t,[r[k] if not r['saturated'] else np.nan for r in rr],color=c)
            sat=[r for r in rr if r['completed'] and r['saturated']]
            a.plot([r['Tm'] for r in sat],[r[k] for r in sat],'x',color=c,ms=6)
    ax[0].set(ylabel='Peak swing (deg)',title='Solid: nonlinear; dashed: linear');ax[0].legend(fontsize=8)
    ax[1].set(ylabel='Absolute peak difference (deg)',title='Crosses: saturated nonlinear runs')
    ax[2].set(ylabel='Relative peak difference (%)',title='Relative to nonlinear peak');ax[2].axhline(5,color='k',ls='--',lw=1)
    for a in ax:a.set(xlabel='Maneuver duration (s)',xlim=(2,15))
    ax[1].set_yscale('log');ax[2].set_yscale('log')
    save(fig,'q3_linear_comparison')

if __name__=='__main__':main()
