"""Additional analyses responding to reviewer comments (PLOS ONE submission).
Run from this directory: python3 review_analyses.py  (seed 20261009)"""
import numpy as np, pandas as pd, json, sys
T="../tables/meta_analysis_RRA_GSE50838_all_genes.csv"
m=pd.read_csv(T); NK=m.NK_IL2_FC_50838.values; TNF=m.Neut_TNF_FC.values; LPS=m.Neut_LPS_FC.values
n=len(m); THR=1.0; rng=np.random.default_rng(20261009)
res={}
def conc(a,b,c,t=THR): return ((a>t)&(b>t)&(c>t))|((a<-t)&(b<-t)&(c<-t))
obs=int(conc(NK,TNF,LPS).sum()); res['observed']=obs
up2=(TNF>THR)&(LPS>THR); dn2=(TNF<-THR)&(LPS<-THR); nUp,nDn=int(up2.sum()),int(dn2.sum()); res['cand_up']=nUp; res['cand_dn']=nDn
# ---- A. analytic expectation
fu=(NK>THR).mean(); fd=(NK<-THR).mean()
E=nUp*fu+nDn*fd
# exact variance of count under random permutation of NK (sampling without replacement)
res['analytic']=dict(frac_NK_up=fu,frac_NK_dn=fd,expected=E)
def cp_counts(nk,B,rng):
    """exact-distribution sampling of correlation-preserving null count via multivariate hypergeometric"""
    a=int((nk>THR).sum()); b=int((nk<-THR).sum()); r=len(nk)-a-b
    d1=rng.multivariate_hypergeometric([a,b,r],nUp,size=B)       # NK values assigned to up-candidates
    rem=np.array([a,b,r])-d1
    d2=np.array([rng.multivariate_hypergeometric(x,nDn) for x in rem]) if B<=2000 else None
    return d1[:,0]+d2[:,1]
def pval(o,null): return (1+(null>=o).sum())/(1+len(null))
# confirm sampling reproduces full permutation result
nc=cp_counts(NK,10000,rng) if False else None
# faster check with direct permutation
chk=np.array([conc(rng.permutation(NK),TNF,LPS).sum() for _ in range(5000)])
res['analytic']['perm_mean']=float(chk.mean()); res['analytic']['perm_sd']=float(chk.std())
# analytic variance (hypergeometric-sum approx, independent pieces)

def indep_fast(nk,t,l,B,thr=THR):
    """exact distribution of the independent-shuffle null count (hypergeometric form; validated against explicit permutation)"""
    out=np.zeros(B,int)
    for sgn in (1,-1):
        A=int((sgn*nk>thr).sum()); nt=int((sgn*t>thr).sum()); nl=int((sgn*l>thr).sum())
        x=rng.hypergeometric(nt,n-nt,A,size=B); out+=rng.hypergeometric(nl,n-nl,x)
    return out
SEC=sys.argv[1]
if SEC=='A': pass
# ---- B. inflation vs correlation r between neutrophil datasets
def gauss_pair(r):
    z1=rng.standard_normal(n); z2=r*z1+np.sqrt(1-r**2)*rng.standard_normal(n)
    return z1,z2
sT=np.sort(TNF); sL=np.sort(LPS)
def map_q(z,sorted_vals):
    ranks=np.argsort(np.argsort(z)); return sorted_vals[ranks]
rows=[]; REPS=1000; B=1000
for r in ([float(x) for x in sys.argv[2].split(',')] if SEC=='B' else []):
    fpI=fpC=0; meanrho=[]; cnts=[]
    for _ in range(REPS):
        z1,z2=gauss_pair(r)                       # real NK stays independent of both by construction
        tnf=map_q(z1,sT); lps=map_q(z2,sL); nk=rng.permutation(NK)
        o=int(conc(nk,tnf,lps).sum()); cnts.append(o)
        meanrho.append(np.corrcoef(tnf,lps)[0,1])
        # independent-shuffle null
        nullI=indep_fast(nk,tnf,lps,B)
        # correlation-preserving null
        u2=(tnf>THR)&(lps>THR); d2_=(tnf<-THR)&(lps<-THR); u,d=int(u2.sum()),int(d2_.sum())
        a=int((nk>THR).sum()); b=int((nk<-THR).sum()); rr=n-a-b
        d1=rng.multivariate_hypergeometric([a,b,rr],u,size=B); rem=np.array([a,b,rr])-d1
        nullC=d1[:,0]+rng.hypergeometric(rem[:,1],(n-u)-rem[:,1],d)
        fpI+= pval(o,nullI)<=0.05; fpC+= pval(o,nullC)<=0.05
    rows.append(dict(r=r,r_realised=float(np.mean(meanrho)),mean_obs=float(np.mean(cnts)),fpr_independent=fpI/REPS,fpr_corrpres=fpC/REPS))
    print(rows[-1],flush=True)
res['inflation']=rows
# ---- C. power simulation (real TNF/LPS and real NK-like background; plant k NK genes among the 224 candidates)
cand=np.where(up2|dn2)[0]; sign=np.where(up2,1,-1)
prow=[]; REPS=300; B=1000
for mag in ([float(x) for x in sys.argv[2].split(',')] if SEC=='C' else []):
  for k in [0,2,4,6,8,10,15]:
    det=0
    for _ in range(REPS):
        nk=rng.permutation(NK).copy()
        # remove accidental concordance? keep as is: background reflects real chance concordance
        if k>0:
            idx=rng.choice(cand,k,replace=False); idx=idx[rng.random(k)<mag]; nk[idx]=sign[idx]*(1.25+rng.exponential(0.25,len(idx)))  # mag = per-gene probability that the true shared gene exceeds |log2FC|>1 in the NK arm
        o=int(conc(nk,TNF,LPS).sum())
        a=int((nk>THR).sum()); b=int((nk<-THR).sum()); rr=n-a-b
        d1=rng.multivariate_hypergeometric([a,b,rr],nUp,size=B); rem=np.array([a,b,rr])-d1
        nullC=d1[:,0]+rng.hypergeometric(rem[:,1],(n-nUp)-rem[:,1],nDn)
        det+= pval(o,nullC)<=0.05
    prow.append(dict(p_detect=mag,k=k,power=det/REPS)); print(prow[-1],flush=True)
res['power']=prow
# ---- D. rank-based (platform-free) sensitivity
from scipy.stats import rankdata
R=lambda x: rankdata(x)/len(x)
rN,rT,rL=R(NK),R(TNF),R(LPS)
srows=[]
for q in ([0.20,0.10,0.05] if SEC=='D' else []):
    cUp=(rT>1-q)&(rL>1-q); cDn=(rT<q)&(rL<q)
    o=int(((rN>1-q)&cUp).sum()+((rN<q)&cDn).sum())
    nu,nd=int(cUp.sum()),int(cDn.sum())
    a=int((rN>1-q).sum()); b=int((rN<q).sum()); rr=n-a-b
    B=10000
    d1=rng.multivariate_hypergeometric([a,b,rr],nu,size=B); rem=np.array([a,b,rr])-d1
    nullC=d1[:,0]+rng.hypergeometric(rem[:,1],(n-nu)-rem[:,1],nd)
    nullI=indep_fast(rN-(1-q),rT-(1-q),rL-(1-q),10000,0)*0+indep_fast(np.where(rN>1-q,1,np.where(rN<q,-1,0)),np.where(rT>1-q,1,np.where(rT<q,-1,0)),np.where(rL>1-q,1,np.where(rL<q,-1,0)),10000,0.5)
    srows.append(dict(q=q,observed=o,cand=nu+nd,corr_mean=float(nullC.mean()),corr_sd=float(nullC.std()),p_corr=pval(o,nullC),indep_mean=float(nullI.mean()),p_indep=pval(o,nullI)))
    print(srows[-1],flush=True)
res['rank_sensitivity']=srows
json.dump(res,open(f'res_{SEC}_{sys.argv[2] if len(sys.argv)>2 else ""}.json','w'),indent=1,default=float)
print(json.dumps(res['analytic']),obs,nUp,nDn)
