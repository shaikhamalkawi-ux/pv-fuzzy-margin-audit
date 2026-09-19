#!/usr/bin/env python3
from pathlib import Path
import json, math, hashlib
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
RES=ROOT/'results'
checks=[]

def ck(name, cond, detail=''):
    checks.append((name,bool(cond),detail))

def close(a,b,tol=1e-9):
    return math.isfinite(float(a)) and abs(float(a)-float(b))<=tol

df=pd.read_csv(RES/'R1_FULL_720_PROFILE_BOUNDARIES.csv')
key=json.loads((RES/'R1_KEY_RESULTS.json').read_text())
alpha=pd.read_csv(RES/'R1_GLOBAL_ALPHA_METRICS.csv')
crisp=pd.read_csv(RES/'R1_CRISP_VS_FUZZY_SUPPORT.csv')
poss=pd.read_csv(RES/'R1_GUARD_POSSIBILITY_BY_WINDOW.csv')
eq=pd.read_csv(RES/'R1_EQUAL_ORDER_FUZZY_FIDELITY.csv')
sens=pd.read_csv(RES/'R1_MEMBERSHIP_POSSIBILITY_SENSITIVITY.csv')

ck('720 profiles',len(df)==720,str(len(df)))
ck('12 windows',df.window_id.nunique()==12,str(df.window_id.nunique()))
for m in ['FOM','SPPK35','SPPK45','CM-SPPK45']:
    ck(f'{m} 720 unique crossings',(df[m+'_crossings']==1).sum()==720,str((df[m+'_crossings']==1).sum()))
ck('Aggregate15 710 unique crossings',(df['Aggregate15_crossings']==1).sum()==710,str((df['Aggregate15_crossings']==1).sum()))
cm=df['CM-SPPK45_required_guard']
ck('CM45 max required guard',close(100*cm.max(),5.57151114230254,1e-9),f'{100*cm.max():.12f}%')
ck('CM45 >2.5 count',(cm>0.025).sum()==9,str((cm>0.025).sum()))
ck('all exceedances in W11',set(df.loc[cm>0.025,'window_id'])=={'W11'},str(set(df.loc[cm>0.025,'window_id'])))
worst=df.loc[cm.idxmax()]
ck('worst timestamp',worst.sample_timestamp=='2010-07-30 13:07:35',str(worst.sample_timestamp))
ck('worst FOM SCR',close(worst.FOM_SCR,0.9859028882450858,1e-10),f'{worst.FOM_SCR:.12f}')
ck('worst CM45 SCR',close(worst['CM-SPPK45_SCR'],0.9338721001314096,1e-10),f"{worst['CM-SPPK45_SCR']:.12f}")
ck('crisp max guard',close(100*crisp.crisp_minute_mean_required_guard.max(),0.8687140424121701,1e-9),f'{100*crisp.crisp_minute_mean_required_guard.max():.12f}%')
w11p=poss[poss.window_id=='W11'].iloc[0]
ck('primary possibility exceedance',close(w11p.possibility_exceed,1.0),str(w11p.possibility_exceed))
ck('primary necessity exceedance',close(w11p.necessity_exceed,0.0),str(w11p.necessity_exceed))
cmalpha=alpha[alpha.model=='CM-SPPK45'].set_index('alpha')
ck('alpha0 max guard',close(100*cmalpha.loc[0.0,'max_guard_req'],5.57151114230254,1e-9),f"{100*cmalpha.loc[0.0,'max_guard_req']:.12f}%")
ck('alpha1 max guard',close(100*cmalpha.loc[1.0,'max_guard_req'],5.532085682900823,1e-8),f"{100*cmalpha.loc[1.0,'max_guard_req']:.12f}%")
ck('equal-order alpha0 gain',close(eq.loc[eq.alpha==0,'CM45_vs_SPPK45_improvement_pct'].iloc[0],36.380178,1e-5),str(eq.loc[eq.alpha==0,'CM45_vs_SPPK45_improvement_pct'].iloc[0]))
sp=dict(zip(sens.variant,sens.possibility_guard_exceedance))
ck('triangular possibility sensitivity',close(sp['TRI50_MIN'],0.493640,1e-6),str(sp['TRI50_MIN']))
ck('depth possibility sensitivity',close(sp['DEPTH_PI'],0.733333,1e-6),str(sp['DEPTH_PI']))

bad=[x for x in checks if not x[1]]
for name,ok,detail in checks: print(('PASS' if ok else 'FAIL'),name,'-',detail)
print(f'\nSUMMARY: {len(checks)-len(bad)}/{len(checks)} checks PASS')
raise SystemExit(1 if bad else 0)
