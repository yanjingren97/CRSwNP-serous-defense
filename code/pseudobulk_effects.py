from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, spearmanr

BASE = ROOT / "results" / "pseudobulk"
CELL_MIN = 50

CONTRASTS = [
    ("GSE235711", "Ethmoid", "CRSsNP", "Healthy"),
    ("GSE235711", "Ethmoid", "CRSwNP", "Healthy"),
    ("GSE276503", "Inferior_turbinate", "CRSwNP", "Healthy"),
    ("GSE261706", "Inferior_turbinate", "AR", "Healthy"),
]


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    x = p[ok]; order = np.argsort(x); ranked = x[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked)+1)
    q = np.minimum.accumulate(q[::-1])[::-1].clip(0,1)
    inv = np.empty_like(order); inv[order] = np.arange(len(order))
    out[np.flatnonzero(ok)] = q[inv]
    return out


def load(ds):
    counts = pd.read_csv(BASE / f"{ds}_pseudobulk_counts.tsv.gz", sep="\t", index_col=0)
    meta = pd.read_csv(BASE / f"{ds}_pseudobulk_metadata.csv").set_index("column").loc[counts.columns]
    return counts, meta


def hedges_rows(a, b):
    n1, n0 = a.shape[1], b.shape[1]
    sp = np.sqrt(((n1-1)*np.var(a,axis=1,ddof=1)+(n0-1)*np.var(b,axis=1,ddof=1))/(n1+n0-2))
    d = (np.mean(a,axis=1)-np.mean(b,axis=1))/np.where(sp==0,np.nan,sp)
    return d*(1-3/(4*(n1+n0)-9))


def main():
    cache = {ds: load(ds) for ds in {x[0] for x in CONTRASTS}}
    results=[]; summaries=[]
    for ds,tissue,case,ctrl in CONTRASTS:
        counts, meta = cache[ds]
        for ct in sorted(meta.celltype.unique()):
            m=meta[(meta.tissue==tissue)&(meta.celltype==ct)&(meta.disease.isin([case,ctrl]))&(meta.n_cells>=CELL_MIN)]
            nc=(m.disease==case).sum(); nn=(m.disease==ctrl).sum()
            if nc<2 or nn<2: continue
            x=counts[m.index]
            lib=x.sum(axis=0).replace(0,np.nan)
            cpm=x.div(lib,axis=1)*1e6
            keep=(cpm>=1).sum(axis=1)>=2
            y=np.log2(cpm.loc[keep]+0.5)
            a=y.loc[:,m.index[m.disease==case]].to_numpy()
            b=y.loc[:,m.index[m.disease==ctrl]].to_numpy()
            stat,p=ttest_ind(a,b,axis=1,equal_var=False,nan_policy="omit")
            tab=pd.DataFrame({"dataset":ds,"tissue":tissue,"case":case,"control":ctrl,"celltype":ct,
                              "gene":y.index,"n_case":nc,"n_control":nn,
                              "log2FC":np.mean(a,axis=1)-np.mean(b,axis=1),
                              "hedges_g":hedges_rows(a,b),"welch_p":p})
            tab["fdr_bh"]=bh(tab.welch_p)
            results.append(tab)
            summaries.append({"dataset":ds,"tissue":tissue,"case":case,"celltype":ct,
                              "n_case":nc,"n_control":nn,"genes_tested":len(tab),
                              "fdr10_absfc05":int(((tab.fdr_bh<.1)&(tab.log2FC.abs()>=.5)).sum()),
                              "p05_absfc05":int(((tab.welch_p<.05)&(tab.log2FC.abs()>=.5)).sum())})
    res=pd.concat(results,ignore_index=True)
    res.to_csv(BASE/"all_pseudobulk_gene_effects.csv.gz",index=False,compression="gzip")
    pd.DataFrame(summaries).to_csv(BASE/"pseudobulk_effect_summary.csv",index=False)

    # Cross-study CRSwNP concordance, without pooling sites or samples.
    concord=[]
    a=res[(res.dataset=="GSE235711")&(res.case=="CRSwNP")]
    b=res[(res.dataset=="GSE276503")&(res.case=="CRSwNP")]
    for ct in sorted(set(a.celltype)&set(b.celltype)):
        z=a[a.celltype==ct][["gene","log2FC"]].merge(
            b[b.celltype==ct][["gene","log2FC"]],on="gene",suffixes=("_ethmoid","_IT"))
        rho,p=spearmanr(z.log2FC_ethmoid,z.log2FC_IT)
        strong=z[(z.log2FC_ethmoid.abs()>=.5)&(z.log2FC_IT.abs()>=.5)]
        concord.append({"celltype":ct,"genes_common":len(z),"spearman_rho":rho,"spearman_p":p,
                        "strong_both":len(strong),"strong_same_direction":int((np.sign(strong.log2FC_ethmoid)==np.sign(strong.log2FC_IT)).sum()),
                        "strong_direction_rate":np.mean(np.sign(strong.log2FC_ethmoid)==np.sign(strong.log2FC_IT)) if len(strong) else np.nan})
    pd.DataFrame(concord).to_csv(BASE/"crswnp_cross_study_concordance.csv",index=False)

    # Candidate genes shared in direction across AR and both CRSwNP contrasts.
    ar=res[(res.dataset=="GSE261706")&(res.case=="AR")]
    shared=[]
    for ct in sorted(set(ar.celltype)&set(a.celltype)&set(b.celltype)):
        z=ar[ar.celltype==ct][["gene","log2FC"]].rename(columns={"log2FC":"AR_log2FC"})
        z=z.merge(a[a.celltype==ct][["gene","log2FC"]].rename(columns={"log2FC":"CRSwNP_ethmoid_log2FC"}),on="gene")
        z=z.merge(b[b.celltype==ct][["gene","log2FC"]].rename(columns={"log2FC":"CRSwNP_IT_log2FC"}),on="gene")
        fc=z[["AR_log2FC","CRSwNP_ethmoid_log2FC","CRSwNP_IT_log2FC"]]
        same=(np.sign(fc).nunique(axis=1)==1)
        z["celltype"]=ct; z["same_direction_all3"]=same
        z["min_abs_log2FC"]=fc.abs().min(axis=1)
        z["mean_abs_log2FC"]=fc.abs().mean(axis=1)
        shared.append(z.sort_values(["same_direction_all3","min_abs_log2FC"],ascending=False))
    pd.concat(shared,ignore_index=True).to_csv(BASE/"ar_crswnp_shared_direction_candidates.csv.gz",index=False,compression="gzip")
    print(pd.DataFrame(summaries).to_string(index=False))
    print("\nCONCORDANCE\n",pd.DataFrame(concord).to_string(index=False))


if __name__ == "__main__": main()
