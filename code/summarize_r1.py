#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np, pandas as pd
from scipy.stats import spearmanr
MODELS=['Aggregate15','SPPK35','SPPK45','CM-SPPK45']
ALPHAS=[0,0.25,0.5,0.75,1.0]

def union_length(intervals):
    x=sorted((a,b) for a,b in intervals if np.isfinite(a) and np.isfinite(b) and b>a)
    if not x:return 0.0
    s,e=x[0]; tot=0.
    for a,b in x[1:]:
        if a<=e:e=max(e,b)
        else:tot+=e-s;s,e=a,b
    return tot+e-s

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--boundaries',required=True)
    ap.add_argument('--crisp-boundaries',required=True)
    ap.add_argument('--window-selection',required=True)
    ap.add_argument('--outdir',required=True)
    a=ap.parse_args();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(a.boundaries); crispb=pd.read_csv(a.crisp_boundaries); win=pd.read_csv(a.window_selection)
    rows=[]
    for m in MODELS:
        v=df[np.isfinite(df.FOM_SCR)&np.isfinite(df[m+'_SCR'])]
        e=100*v[m+'_margin_rel_error']; gr=100*v[m+'_required_guard']
        rows.append(dict(model=m,n=len(v),mean_error_pct=e.mean(),median_error_pct=e.median(),p95_error_pct=e.quantile(.95),max_error_pct=e.max(),
                         false_stable=int((v[m+'_SCR']<v.FOM_SCR).sum()),false_unstable=int((v[m+'_SCR']>v.FOM_SCR).sum()),
                         max_required_guard_pct=gr.max(),p95_required_guard_pct=gr.quantile(.95),guard_exceed_2p5=int((gr>2.5).sum())))
    pd.DataFrame(rows).to_csv(out/'R1_PROFILE_MODEL_SUMMARY_REPLAY.csv',index=False)
    rows=[]
    for wid,g in df.groupby('window_id',sort=False):
        for al in ALPHAS:
            sub=g[g.membership>=al-1e-12]
            for m in MODELS:
                v=sub[np.isfinite(sub.FOM_SCR)&np.isfinite(sub[m+'_SCR'])]
                if not len(v):continue
                fl,fh=v.FOM_SCR.min(),v.FOM_SCR.max(); rl,rh=v[m+'_SCR'].min(),v[m+'_SCR'].max(); fm=.5*(fl+fh)
                inter=max(0,min(fh,rh)-max(fl,rl));uni=max(fh,rh)-min(fl,rl)
                fs=[(r,f) for r,f in zip(v[m+'_SCR'],v.FOM_SCR) if r<f]
                fu=[(f,r) for r,f in zip(v[m+'_SCR'],v.FOM_SCR) if r>f]
                rows.append(dict(window_id=wid,alpha=al,n=len(v),model=m,FOM_lo=fl,FOM_hi=fh,ROM_lo=rl,ROM_hi=rh,
                                 hausdorff_rel=max(abs(rl-fl),abs(rh-fh))/max(abs(fm),1e-12),overlap=inter/uni if uni else 1.,
                                 fs_profiles=len(fs),fu_profiles=len(fu),fs_union_width=union_length(fs),fu_union_width=union_length(fu),
                                 guard_req_lo=v[m+'_required_guard'].min(),guard_req_hi=v[m+'_required_guard'].max()))
    ac=pd.DataFrame(rows); ac.to_csv(out/'R1_ALPHA_CUT_SUMMARY_REPLAY.csv',index=False)
    gr=ac.groupby(['model','alpha']).agg(windows=('window_id','nunique'),mean_hausdorff_rel=('hausdorff_rel','mean'),median_hausdorff_rel=('hausdorff_rel','median'),
                                        max_hausdorff_rel=('hausdorff_rel','max'),mean_overlap=('overlap','mean'),min_overlap=('overlap','min'),max_guard_req=('guard_req_hi','max'),
                                        mean_guard_req_hi=('guard_req_hi','mean'),total_fs_profiles=('fs_profiles','sum'),total_fu_profiles=('fu_profiles','sum'),
                                        mean_fs_union_width=('fs_union_width','mean'),max_fs_union_width=('fs_union_width','max')).reset_index()
    gr.to_csv(out/'R1_GLOBAL_ALPHA_METRICS_REPLAY.csv',index=False)
    prow=[]
    for wid,g in df.groupby('window_id',sort=False):
        ex=g['CM-SPPK45_required_guard']>.025; safe=~ex
        pe=float(g.loc[ex,'membership'].max()) if ex.any() else 0.0
        ps=float(g.loc[safe,'membership'].max()) if safe.any() else 0.0
        prow.append(dict(window_id=wid,n_exceed=int(ex.sum()),max_required_guard_pct=100*g['CM-SPPK45_required_guard'].max(),
                         possibility_exceed=pe,necessity_exceed=1-ps,possibility_safe=ps,necessity_safe=1-pe))
    pd.DataFrame(prow).to_csv(out/'R1_GUARD_POSSIBILITY_BY_WINDOW_REPLAY.csv',index=False)
    support=df.groupby('window_id')['CM-SPPK45_required_guard'].max(); c=crispb.set_index('window_id')['CM-SPPK45_required_guard']
    cv=win.set_index('window_id')
    comp=pd.DataFrame({'window_id':support.index,'crisp_minute_mean_required_guard':c.reindex(support.index).values,'fuzzy_support_max_required_guard':support.values})
    for col in ['cv_band','mean_bin','mean_GHI','spatial_cv','temporal_mean_cv']:
        comp[col]=cv.reindex(support.index)[col].values
    comp['understatement_pp']=100*(comp.fuzzy_support_max_required_guard-comp.crisp_minute_mean_required_guard)
    comp.to_csv(out/'R1_CRISP_VS_FUZZY_SUPPORT_REPLAY.csv',index=False)
    ar=[]
    for p in ['spatial_cv','temporal_mean_cv','mean_GHI']:
        rho,pv=spearmanr(comp[p],comp.fuzzy_support_max_required_guard)
        ar.append(dict(predictor=p,spearman_rho=rho,p_value=pv,n=len(comp)))
    pd.DataFrame(ar).to_csv(out/'R1_WINDOW_VARIABILITY_ASSOCIATION_REPLAY.csv',index=False)
    print('R1 summaries regenerated from locked boundary files.')
if __name__=='__main__':main()
