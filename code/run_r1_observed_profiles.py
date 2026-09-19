#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import argparse, json, sys, os
import numpy as np, pandas as pd
from scipy.optimize import brentq
from concurrent.futures import ProcessPoolExecutor, as_completed

MODELS=['Aggregate15','SPPK35','SPPK45','CM-SPPK45']
ORDER=['DH1','DH9','DH8','DH2','DH4','DH3','AP1','AP3','AP4','AP6']

_G={}
def init_worker(upstream_root):
    root=Path(upstream_root)
    sys.path.insert(0,str(root/'code'))
    sys.path.insert(0,str(root/'legacy'/'code'))
    from extension_core import equilibrium_family, jacobian, abscissa, reduced, aggregate_basis
    from run_trackB_heterogeneity import string_transverse_projection
    from p001_fast_projected import projected_jacobian_directional
    pk=np.load(root/'inputs/trackB_pk60_basis.npz'); cm=np.load(root/'inputs/trackB_cm_sppk45_profileanchors_basis.npz')
    va,wa=aggregate_basis()
    _G.update(dict(equilibrium_family=equilibrium_family,jacobian=jacobian,abscissa=abscissa,reduced=reduced,
                   string_transverse_projection=string_transverse_projection,projected_jacobian_directional=projected_jacobian_directional,
                   bases={'Aggregate15':(va,wa),'SPPK35':(pk['V'][:,:35],pk['Usp'][:,:35]),
                          'SPPK45':(pk['V'][:,:45],pk['Usp'][:,:45]),'CM-SPPK45':(cm['V'],cm['Usp'])}))

def eval_profile(rec):
    vec=np.array([float(rec[n]) for n in ORDER],float)
    S=np.tile(vec,(40,1))
    ef=_G['equilibrium_family'](S)
    c,y,z,dc=ef(.5)
    Vt,Wt=_G['string_transverse_projection'](c)
    At=_G['projected_jacobian_directional'](y,c,Vt,Wt,rel=1.5e-7)
    at=_G['abscissa'](At)
    cache={}
    def vals(x):
        x=float(x)
        if x in cache:return cache[x]
        c,y,z,dc=ef(x); A=_G['jacobian'](z,c); ac=_G['abscissa'](A)
        o={'FOM':max(float(ac.real),float(at.real)),'fom_freq':abs(float(ac.imag))/(2*np.pi),
           'fom_coherent':float(ac.real),'fom_transverse':float(at.real)}
        for n,(V,U) in _G['bases'].items():o[n]=float(_G['abscissa'](_G['reduced'](A,V,U)).real)
        cache[x]=o; return o
    scan=np.linspace(.15,1.8,12)
    for x in scan:vals(x)
    roots={}; status={}
    for model in ['FOM']+MODELS:
        vv=np.array([cache[float(x)][model] for x in scan])
        ix=np.where(vv[:-1]*vv[1:]<0)[0]
        good=[j for j in ix if vv[j]<0 and vv[j+1]>0]
        status[model+'_crossings']=len(good)
        if len(good)!=1: roots[model]=None
        else:
            j=good[0]
            roots[model]=brentq(lambda xx:vals(xx)[model],float(scan[j]),float(scan[j+1]),xtol=5e-6,rtol=1e-8,maxiter=30)
    out={k:rec[k] for k in ['window_id','window_timestamp','sample_second','sample_timestamp','membership','mean_GHI','spatial_cv']}
    out.update(status)
    for model,r in roots.items():
        out[model+'_Xg']=r
        out[model+'_SCR']=1/r if r else np.nan
    sf=out.get('FOM_SCR',np.nan)
    if np.isfinite(sf):
        for model in MODELS:
            sr=out.get(model+'_SCR',np.nan)
            if np.isfinite(sr):
                out[model+'_margin_rel_error']=abs(sr/sf-1)
                out[model+'_signed_rel_shift']=sr/sf-1
                out[model+'_required_guard']=max(0.0,sf/sr-1)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--tasks',required=True); ap.add_argument('--upstream-root',required=True); ap.add_argument('--out',required=True)
    ap.add_argument('--workers',type=int,default=5)
    a=ap.parse_args()
    df=pd.read_csv(a.tasks); recs=df.to_dict('records'); rows=[None]*len(recs)
    with ProcessPoolExecutor(max_workers=a.workers,initializer=init_worker,initargs=(a.upstream_root,)) as ex:
        futs={ex.submit(eval_profile,r):i for i,r in enumerate(recs)}
        done=0
        for fut in as_completed(futs):
            i=futs[fut]; rows[i]=fut.result(); done+=1
            if done%50==0 or done==len(recs):print(f'{done}/{len(recs)}',flush=True)
    out=pd.DataFrame(rows); Path(a.out).parent.mkdir(parents=True,exist_ok=True); out.to_csv(a.out,index=False)
    print(out[['FOM_crossings']+[m+'_crossings' for m in MODELS]].apply(pd.Series.value_counts).fillna(0).to_string())
if __name__=='__main__':main()
