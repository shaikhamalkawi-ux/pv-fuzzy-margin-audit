#!/usr/bin/env python3
from pathlib import Path
import math
import pandas as pd

ROOT=Path(__file__).resolve().parent.parent
RES=ROOT/'results'
checks=[]

def ck(name,cond,detail=''):
    checks.append((name,bool(cond),detail))

def close(a,b,tol=1e-9):
    return math.isfinite(float(a)) and abs(float(a)-float(b))<=tol

crisp=pd.read_csv(RES/'R1_CRISP_MEAN_BOUNDARIES.csv')
alpha=pd.read_csv(RES/'R1_GLOBAL_ALPHA_METRICS.csv')
poss=pd.read_csv(RES/'R1_GUARD_POSSIBILITY_BY_WINDOW.csv')
base=pd.read_csv(RES/'R1_BASELINE_COMPARISON_SUMMARY.csv')

ck('12 crisp windows',len(crisp)==12,str(len(crisp)))
ck('crisp max guard',close(100*crisp['CM-SPPK45_required_guard'].max(),0.8687140425906792,1e-8),f"{100*crisp['CM-SPPK45_required_guard'].max():.12f}%")
w11=poss[poss.window_id=='W11'].iloc[0]
ck('W11 exceedance count',int(w11.n_exceed)==9,str(int(w11.n_exceed)))
ck('W11 support max guard',close(float(w11.max_required_guard_pct),5.57151114230254,1e-9),f"{w11.max_required_guard_pct:.12f}%")
ck('W11 possibility exceedance',close(w11.possibility_exceed,1.0),str(w11.possibility_exceed))
cm=alpha[alpha.model=='CM-SPPK45'].set_index('alpha')
ck('alpha0 support max',close(100*cm.loc[0.0,'max_guard_req'],5.57151114230254,1e-9),f"{100*cm.loc[0.0,'max_guard_req']:.12f}%")
ck('alpha1 max',close(100*cm.loc[1.0,'max_guard_req'],5.53208568290097,1e-8),f"{100*cm.loc[1.0,'max_guard_req']:.12f}%")
obs=base[base.baseline=='Observed support'].iloc[0]
ck('observed support summary',close(100*obs.max_guard_summary,5.57151114230254,1e-9),f"{100*obs.max_guard_summary:.12f}%")

bad=[x for x in checks if not x[1]]
for name,ok,detail in checks:
    print(('PASS' if ok else 'FAIL'),name,'-',detail)
print(f"SUMMARY: {len(checks)-len(bad)}/{len(checks)} checks PASS")
raise SystemExit(1 if bad else 0)
