import numpy as np, pandas as pd, json
from scipy.stats import rankdata
m=pd.read_csv("../tables/meta_analysis_RRA_GSE50838_all_genes.csv")
NK=m.NK_IL2_FC_50838.values;T=m.Neut_TNF_FC.values;L=m.Neut_LPS_FC.values;n=len(m);rng=np.random.default_rng(20261009)
R=lambda x: rankdata(x)/len(x)
rN,rT,rL=R(NK),R(T),R(L); B=10000; out=[]
pv=lambda o,nl:(1+(nl>=o).sum())/(1+len(nl))
for q in [0.20,0.10,0.05]:
    cU=(rT>1-q)&(rL>1-q); cD=(rT<q)&(rL<q); nu,nd=int(cU.sum()),int(cD.sum())
    o=int(((rN>1-q)&cU).sum()+((rN<q)&cD).sum())
    a=int((rN>1-q).sum()); b=int((rN<q).sum())
    # correlation-preserving null: assign NK labels at random to the nu up-candidates and nd down-candidates
    nullC=np.empty(B,int)
    for i in range(B):
        idx=rng.permutation(n)[:nu+nd]
        nullC[i]=(rN[idx[:nu]]>1-q).sum()+(rN[idx[nu:]]<q).sum()
    # independent null (exact hypergeometric form)
    nullI=np.zeros(B,int)
    for sg in (1,-1):
        A=a if sg==1 else b; nt=int((rT>1-q).sum()) if sg==1 else int((rT<q).sum()); nl=int((rL>1-q).sum()) if sg==1 else int((rL<q).sum())
        x=rng.hypergeometric(nt,n-nt,A,size=B); nullI+=rng.hypergeometric(nl,n-nl,x)
    out.append(dict(q=q,observed=o,candidates=nu+nd,corr_mean=float(nullC.mean()),corr_sd=float(nullC.std()),p_corr=pv(o,nullC),indep_mean=float(nullI.mean()),indep_sd=float(nullI.std()),p_indep=pv(o,nullI)))
    print(out[-1],flush=True)
json.dump(out,open('res_D.json','w'),indent=1,default=float)
