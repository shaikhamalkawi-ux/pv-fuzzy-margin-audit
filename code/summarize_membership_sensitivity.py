#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd, numpy as np
ALPHAS=[0.25,0.5,0.75,1.0]
VARIANTS=['T10_MIN','T10_PRODUCT','T10R_MIN','T25_MIN','TRI50_MIN','DEPTH_PI']

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--boundaries',required=True); ap.add_argument('--membership-values',required=True); ap.add_argument('--outdir',required=True)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    b=pd.read_csv(a.boundaries); m=pd.read_csv(a.membership_values)
    d=b.merge(m,on=['window_id','sample_second'],how='left',validate='one_to_one')
    rows=[]
    for v in VARIANTS:
        col=v+'_NORM'
        for al in ALPHAS:
            for wid,g in d.groupby('window_id',sort=False):
                q=g[g[col]>=al-1e-12]
                if not len(q):continue
                rows.append(dict(variant=v,alpha=al,window_id=wid,n=len(q),max_required_guard=q['CM-SPPK45_required_guard'].max(),max_error=q['CM-SPPK45_margin_rel_error'].max(),exceeds_guard=bool((q['CM-SPPK45_required_guard']>.025).any())))
    detail=pd.DataFrame(rows); detail.to_csv(out/'R1_MEMBERSHIP_SENSITIVITY_BY_WINDOW_REPLAY.csv',index=False)
    summ=detail.groupby(['variant','alpha']).agg(windows_nonempty=('window_id','nunique'),min_n=('n','min'),max_n=('n','max'),max_required_guard=('max_required_guard','max'),mean_window_max_guard=('max_required_guard','mean'),max_error=('max_error','max'),windows_exceeding_guard=('exceeds_guard','sum')).reset_index()
    summ.to_csv(out/'R1_MEMBERSHIP_SENSITIVITY_SUMMARY_NORMALIZED_REPLAY.csv',index=False)
    poss=[]
    for v in VARIANTS:
        col=v+'_NORM'; ex=d['CM-SPPK45_required_guard']>.025; safe=~ex
        pe=float(d.loc[ex,col].max()) if ex.any() else 0.; ps=float(d.loc[safe,col].max()) if safe.any() else 0.
        poss.append(dict(variant=v,possibility_guard_exceedance=pe,necessity_guard_exceedance=1-ps,possibility_guard_sufficiency=ps,necessity_guard_sufficiency=1-pe))
    pd.DataFrame(poss).to_csv(out/'R1_MEMBERSHIP_POSSIBILITY_SENSITIVITY_REPLAY.csv',index=False)
    print('Membership-sensitivity summaries regenerated from normalized archived membership values.')
if __name__=='__main__':main()
