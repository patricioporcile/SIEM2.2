import pandas as pd, numpy as np, json
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
U="../data/GSE288312_Gene_counts_Wright_2025.csv"  # download from GEO (GSE288312 supplementary file)
T="../tables/"
c=pd.read_csv(U,index_col=0).set_index("GeneID")
print(c.shape, c.index.duplicated().sum())
c=c.groupby(level=0).sum()
cols=[x for x in c.columns if x.startswith("UT")]+[x for x in c.columns if x.startswith("IL1B")]
X=c[cols].T
meta=pd.DataFrame({"cond":["UT"]*6+["IL1B"]*3},index=cols)
X=X.loc[:,X.sum()>=10]
dds=DeseqDataSet(counts=X,metadata=meta,design="~cond",refit_cooks=True,quiet=True); dds.deseq2()
st=DeseqStats(dds,contrast=["cond","IL1B","UT"],quiet=True); st.summary()
r=st.results_df; r.to_csv("GSE288312_IL1B_DESeq2.csv"); print(r.shape)
cand=pd.read_csv(T+"candidates_224_TNF_LPS_concordant.csv"); print(cand.columns.tolist(), len(cand)); print(cand.head(3))
T="../tables/"
r=pd.read_csv("GSE288312_IL1B_DESeq2.csv",index_col=0)
cand=pd.read_csv(T+"candidates_224_TNF_LPS_concordant.csv").set_index("gene")
bg=pd.read_csv(T+"meta_analysis_RRA_GSE50838_all_genes.csv").Gene
rng=np.random.default_rng(20261030); B=10000; out={}
def test(pool,label,extra=None):
    pv=r.loc[r.index.isin(pool)]
    lfc=pv.log2FoldChange.dropna()
    cc=cand[cand.index.isin(lfc.index)]
    nu=int((cc.direction=="up").sum()); nd=int((cc.direction=="down").sum())
    f=lfc.loc[cc.index].values; d=np.where(cc.direction=="up",1,-1)
    obs=int(((np.abs(f)>1)&(np.sign(f)==d)).sum())
    n=len(lfc); a=int((lfc>1).sum()); b=int((lfc<-1).sum()); rest=n-a-b
    d1=rng.multivariate_hypergeometric([a,b,rest],nu,size=B); rem=np.array([a,b,rest])-d1
    null=d1[:,0]+rng.hypergeometric(rem[:,1],(n-nu)-rem[:,1],nd)
    p=(1+(null>=obs).sum())/(B+1)
    res=dict(background=label,n_pool=n,n_cand=nu+nd,observed=obs,null_mean=round(float(null.mean()),2),null_sd=round(float(null.std()),2),null_max=int(null.max()),p=float(p))
    print(res); out[label]=res
    return cc[(np.abs(f)>1)&(np.sign(f)==d)].index.tolist()
g1=test(r.index,"unrestricted (expressed genes)")
g2=test(set(bg),"restricted (8,581 shared)")
print(g1)
# FDR-based variant
sig=r[(r.padj<0.05)]
cc=cand[cand.index.isin(sig.index)]
f=sig.loc[cc.index,"log2FoldChange"].values; d=np.where(cc.direction=="up",1,-1)
print("concordant & padj<0.05 & |lfc|>1:", int(((np.abs(f)>1)&(np.sign(f)==d)).sum()), "| padj<.05 any lfc:", int((np.sign(f)==d).sum()))
print("candidates present:",int(cand.index.isin(r.index).sum()))
json.dump(out,open("pc_deseq2.json","w"),indent=1)
