#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import argparse, zipfile, hashlib
import numpy as np, pandas as pd

ORDER=['DH1','DH9','DH8','DH2','DH4','DH3','AP1','AP3','AP4','AP6']
RAW_ORDER=['DH3','DH4','DH5','DH10','DH11','DH9','DH2','DH1','DH1T','AP6','AP6T','AP1','AP3','AP5','AP4','AP7','DH6','DH7','DH8']
MEAN_EDGES=np.array([400.,575.,750.,925.,1100.])

def trap_mu(x,a,b,c,d):
    x=np.asarray(x,float); out=np.zeros_like(x)
    if b>a:
        m=(x>a)&(x<b); out[m]=(x[m]-a)/(b-a)
    out[(x>=b)&(x<=c)]=1.
    if d>c:
        m=(x>c)&(x<d); out[m]=(d-x[m])/(d-c)
    if a==b: out[x==a]=1.
    if c==d: out[x==d]=1.
    return np.clip(out,0,1)

def read_daily(z,member,usecols,names):
    with z.open(member) as f:return pd.read_csv(f,header=None,usecols=usecols,names=names)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('raw_zip'); ap.add_argument('outdir'); ap.add_argument('--pilot-per-window',type=int,default=20)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    raw=Path(a.raw_zip); raw_idx={n:i for i,n in enumerate(RAW_ORDER)}
    usecols=[0,1,2,3]+[4+raw_idx[n] for n in ORDER]; names=['seconds','year','doy','hst']+ORDER
    records=[]; minute_store={}
    with zipfile.ZipFile(raw) as z:
        for member in sorted(n for n in z.namelist() if n.lower().endswith('.txt')):
            df=read_daily(z,member,usecols,names)
            for idx,g in df.groupby(['year','doy','hst'],sort=False):
                if len(g)!=60: continue
                arr=g[ORDER].to_numpy(float)
                if not np.isfinite(arr).all() or (arr<100).any() or (arr>1400).any():continue
                vec=arr.mean(axis=0); mu=float(vec.mean()); sd=float(vec.std(ddof=0))
                if not (400<=mu<=1100) or sd<=0:continue
                year,doy,hst=idx; h=int(hst)//100; m=int(hst)%100
                ts=pd.to_datetime(f'{int(year)}-{int(doy):03d}',format='%Y-%j')+pd.Timedelta(hours=h,minutes=m)
                records.append({'timestamp':ts,'mean_GHI':mu,'spatial_cv':sd/mu,
                                'temporal_mean_cv':float(np.mean(arr.std(axis=0,ddof=0)/(arr.mean(axis=0)+1e-12)))})
                minute_store[ts]=(g['seconds'].to_numpy(int),arr)
    cand=pd.DataFrame(records).sort_values('timestamp').reset_index(drop=True)
    q1,q2=cand.spatial_cv.quantile([1/3,2/3]).tolist(); cv_edges=[cand.spatial_cv.min()-1e-12,q1,q2,cand.spatial_cv.max()+1e-12]
    selected=[]
    for ci in range(3):
        lo,hi=cv_edges[ci],cv_edges[ci+1]; subcv=cand[(cand.spatial_cv>=lo)&(cand.spatial_cv<hi)]; cv_center=float(subcv.spatial_cv.median())
        for mi in range(4):
            ml,mh=MEAN_EDGES[mi],MEAN_EDGES[mi+1]; sub=subcv[(subcv.mean_GHI>=ml)&(subcv.mean_GHI<(mh if mi<3 else mh+1e-9))].copy()
            mcenter=(ml+mh)/2
            sub['selection_score']=((sub.mean_GHI-mcenter)/(mh-ml))**2+((sub.spatial_cv-cv_center)/max(hi-lo,1e-12))**2
            r=sub.sort_values(['selection_score','timestamp']).iloc[0]
            selected.append({'window_id':f'W{len(selected)+1:02d}','cv_band':ci+1,'mean_bin':mi+1,'timestamp':r.timestamp,
                             'mean_bin_low':ml,'mean_bin_high':mh,'cv_band_low':lo,'cv_band_high':hi,'mean_GHI':r.mean_GHI,
                             'spatial_cv':r.spatial_cv,'temporal_mean_cv':r.temporal_mean_cv,'selection_score':r.selection_score,'eligible_cell_n':len(sub)})
    sel=pd.DataFrame(selected); sel['timestamp']=sel.timestamp.dt.strftime('%Y-%m-%d %H:%M:%S'); sel.to_csv(out/'R0_WINDOW_SELECTION.csv',index=False)
    tasks=[]; pars=[]
    for r in sel.itertuples(index=False):
        ts=pd.Timestamp(r.timestamp); secs,arr=minute_store[ts]; mus=[]; p={'window_id':r.window_id,'timestamp':r.timestamp}
        for j,n in enumerate(ORDER):
            x=arr[:,j]; aa=float(x.min()); b=float(np.quantile(x,.10)); c=float(np.quantile(x,.90)); d=float(x.max())
            p.update({f'{n}_a_min':aa,f'{n}_b_q10':b,f'{n}_c_q90':c,f'{n}_d_max':d}); mus.append(trap_mu(x,aa,b,c,d))
        pars.append(p); joint=np.min(np.column_stack(mus),axis=1)
        for sec,vec,mu in zip(secs,arr,joint):
            rec={'window_id':r.window_id,'window_timestamp':r.timestamp,'sample_second':int(sec),
                 'sample_timestamp':(ts+pd.Timedelta(seconds=int(sec))).strftime('%Y-%m-%d %H:%M:%S'),
                 'membership':float(mu),'mean_GHI':float(vec.mean()),'spatial_cv':float(vec.std(ddof=0)/vec.mean())}
            rec.update({n:float(v) for n,v in zip(ORDER,vec)}); tasks.append(rec)
    full=pd.DataFrame(tasks); pd.DataFrame(pars).to_csv(out/'R0_FUZZY_INPUT_PARAMETERS.csv',index=False); full.to_csv(out/'R0_OBSERVED_PROFILE_TASKS.csv',index=False)
    pilot=[]
    for wid,g in full.groupby('window_id',sort=False):
        gg=g.sort_values(['membership','sample_second'],ascending=[True,True]).reset_index(drop=True)
        idx=np.unique(np.round(np.linspace(0,len(gg)-1,a.pilot_per_window)).astype(int)); e=gg.iloc[idx].copy(); e['pilot_rank_index']=idx; pilot.append(e)
    pd.concat(pilot,ignore_index=True).to_csv(out/'R0_PILOT_PROFILE_TASKS_240.csv',index=False)
    for name in ['R0_WINDOW_SELECTION.csv','R0_FUZZY_INPUT_PARAMETERS.csv','R0_OBSERVED_PROFILE_TASKS.csv','R0_PILOT_PROFILE_TASKS_240.csv']:
        p=out/name; (out/(name+'.sha256')).write_text(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {name}\n')
    print('eligible minutes',len(cand),'windows',len(sel),'tasks',len(full),'pilot',len(pd.concat(pilot)))
if __name__=='__main__':main()
