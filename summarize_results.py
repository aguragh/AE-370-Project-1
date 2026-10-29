"""Print numerical summaries used in the results text; no simulation rerun."""
import json
import numpy as np
from results_figures import read,group,ROOT

def main():
 r=read('sweep.csv');b=read('boundary_checks.csv');lin=read('linear_comparison.csv')
 tr=json.loads((ROOT/'data'/'transitions.json').read_text())
 print('N runs',len(r),'completed',sum(x['completed'] for x in r),'acceptable',sum(x['acceptable'] for x in r),'stopped',sum(not x['completed'] for x in r))
 print('Transitions',len(tr),'reverse',sum(x['left_acceptable'] for x in tr))
 for m in [.2,.5,1.]:
  for l in [.5,1.75,3.]:
   tt=[t for t in tr if abs(t['m']-m)<1e-8 and t['l']==l]
   print('BOUNDARY',m,l,tt)
 print('BOUNDARY peak max change deg',max(np.rad2deg(abs(x['phi_peak']-y['phi_peak'])) for x,y in zip(b[::2],b[1::2])))
 print('BOUNDARY classification differences',sum(x['acceptable']!=y['acceptable'] for x,y in zip(b[::2],b[1::2])))
 acc=[x for x in b if x['acceptable']]
 for k in ['Tc_min','z_max','pos_err','speed','phi_peak']:print('BOUNDARY accepted metric',k,'range',min(x[k] for x in acc),max(x[k] for x in acc))
 print('BOUNDARY saturation',sum(x['saturated'] for x in acc))
 for l in [.5,1.75,3.]:
  rr=group(r,.5,l)
  for T in [7.,7.5,8.,8.25,8.5,8.75,9.,9.5,10.,12.,15.]:
   x=next(x for x in rr if x['Tm']==T)
   print('NOMINAL',l,T,'peak',np.rad2deg(x['phi_peak']),'res',np.rad2deg(x['phi_res']),'ok',x['acceptable'])
  ll=[x for x in lin if x['l']==l]
  print('LINEAR',l,'unsat >5 times',[x['Tm'] for x in ll if x['completed'] and not x['saturated'] and x['relative_percent']>5])
  print('LINEAR',l,'sat times',[x['Tm'] for x in ll if x['completed'] and x['saturated']])
  for T in [4.,6.,8.,10.,15.]:print('LINEAR SAMPLE',next(x for x in ll if x['Tm']==T))
 corners=read('corner_checks.csv')
 changes=[np.rad2deg(abs(a['phi_peak']-b['phi_peak'])) for a,b in zip(corners[::2],corners[1::2]) if a['completed'] and b['completed']]
 print('CORNERS max change',max(changes),'class mismatches',sum(a['acceptable']!=b['acceptable'] for a,b in zip(corners[::2],corners[1::2])))
 windows=read('window_checks.csv')
 for a,b in zip(windows[::2],windows[1::2]):print('WINDOW',a['l'],a['Tm'],'peak/res diff deg',np.rad2deg(a['phi_peak']-b['phi_peak']),np.rad2deg(a['phi_res']-b['phi_res']),'ok',a['acceptable'],b['acceptable'])
 print('SWEEP runtime median',np.median([x['runtime_s'] for x in r]),'sum',sum(x['runtime_s'] for x in r))

if __name__=='__main__':main()
