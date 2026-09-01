from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import t as tdist, norm, chi2
from secretory_subtype_pseudobulk_validation import load, DEFINITION_GENES

OUT=ROOT/"results/frozen_v1"; OUT.mkdir(parents=True,exist_ok=True)

def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(np.nan_to_num(p,nan=1.)); q=np.full(n,np.nan); last=1.
    for rank,i in reversed(list(enumerate(o,1))):
        if np.isnan(p[i]): continue
        last=min(last,p[i]*n/rank); q[i]=last
    return q

def paired_model(x,m):
    z=m[(m.disease=="CRSwNP")&(m.subtype=="Serous_glandular")&m.tissue.isin(["Nasal_polyp","Ethmoid"])&(m.n_cells>=20)].copy()
    z["column_id"]=z.index
    w=z.reset_index(drop=True).pivot(index="donor",columns="tissue",values="column_id").dropna()
    a,b=w.Nasal_polyp.tolist(),w.Ethmoid.tolist()
    y=x[~x.index.isin(DEFINITION_GENES)]
    keep=((y[a]>0).sum(axis=1)+(y[b]>0).sum(axis=1))>=2; y=y.loc[keep]
    d=y[a].to_numpy()-y[b].to_numpy(); n=d.shape[1]
    eff=d.mean(axis=1); sd=d.std(axis=1,ddof=1); se=sd/np.sqrt(n); tv=np.divide(eff,se,out=np.zeros_like(eff),where=se>0)
    p=2*tdist.sf(np.abs(tv),n-1)
    return pd.DataFrame({"gene":y.index,"effect_paired":eff,"se_paired":se,"ci_low_paired":eff-tdist.ppf(.975,n-1)*se,
                         "ci_high_paired":eff+tdist.ppf(.975,n-1)*se,"p_paired":p,"fdr_paired":bh(p),"n_paired":n})

def unpaired_model(x,m):
    z=m[(m.disease=="CRSwNP")&(m.subtype=="Serous_glandular")&m.tissue.isin(["Nasal_polyp","Inferior_turbinate"])&(m.n_cells>=20)]
    a=z[z.tissue=="Nasal_polyp"].index.tolist(); b=z[z.tissue=="Inferior_turbinate"].index.tolist()
    y=x[~x.index.isin(DEFINITION_GENES)]; keep=((y[a]>0).sum(axis=1)+(y[b]>0).sum(axis=1))>=2; y=y.loc[keep]
    ya,yb=y[a].to_numpy(),y[b].to_numpy(); n1,n0=len(a),len(b)
    eff=ya.mean(axis=1)-yb.mean(axis=1); v1=ya.var(axis=1,ddof=1)/n1; v0=yb.var(axis=1,ddof=1)/n0; se=np.sqrt(v1+v0)
    df=(v1+v0)**2/(v1**2/(n1-1)+v0**2/(n0-1)); tv=np.divide(eff,se,out=np.zeros_like(eff),where=se>0)
    p=2*tdist.sf(np.abs(tv),df); crit=tdist.ppf(.975,df)
    return pd.DataFrame({"gene":y.index,"effect_external":eff,"se_external":se,"ci_low_external":eff-crit*se,
                         "ci_high_external":eff+crit*se,"p_external":p,"fdr_external":bh(p),"n_external_np":n1,"n_external_it":n0})

x1,m1=load("GSE235711"); x2,m2=load("GSE276503")
a=paired_model(x1,m1); b=unpaired_model(x2,m2); z=a.merge(b,on="gene")

# Fixed-effect synthesis is a compact reproducibility summary, not a population-wide clinical estimate.
w1=1/z.se_paired.pow(2); w2=1/z.se_external.pow(2); denom=w1+w2
z["meta_effect"]=(w1*z.effect_paired+w2*z.effect_external)/denom
z["meta_se"]=np.sqrt(1/denom); z["meta_z"]=z.meta_effect/z.meta_se; z["meta_p"]=2*norm.sf(z.meta_z.abs()); z["meta_fdr"]=bh(z.meta_p)
z["meta_ci_low"]=z.meta_effect-1.96*z.meta_se; z["meta_ci_high"]=z.meta_effect+1.96*z.meta_se
z["direction_concordant"]=np.sign(z.effect_paired)==np.sign(z.effect_external)
z["Q"] = w1*(z.effect_paired-z.meta_effect)**2+w2*(z.effect_external-z.meta_effect)**2
z["heterogeneity_p"]=chi2.sf(z.Q,1); z["I2_percent"]=np.maximum(0,(z.Q-1)/z.Q)*100
z.replace([np.inf,-np.inf],np.nan,inplace=True)
z.to_csv(OUT/"formal_serous_gene_models.csv.gz",index=False,compression="gzip")

locked=pd.read_csv(ROOT/"results/locked_40_gene_module.csv").gene
summary=pd.DataFrame([{
 "genes_common":len(z),"direction_concordant_all":z.direction_concordant.mean(),
 "genes_both_effect_le_minus_0_5":int(((z.effect_paired<=-.5)&(z.effect_external<=-.5)).sum()),
 "genes_meta_fdr_0_05":int((z.meta_fdr<.05).sum()),
 "genes_meta_fdr_0_05_concordant_down":int(((z.meta_fdr<.05)&z.direction_concordant&(z.meta_effect<0)).sum()),
 "locked_genes_present":int(z.gene.isin(locked).sum()),
 "locked_genes_concordant_down":int((z[z.gene.isin(locked)].effect_paired.lt(0)&z[z.gene.isin(locked)].effect_external.lt(0)).sum()),
 "locked_median_meta_effect":z[z.gene.isin(locked)].meta_effect.median(),
 "locked_median_I2":z[z.gene.isin(locked)].I2_percent.median()
}])
summary.to_csv(OUT/"formal_serous_gene_model_summary.csv",index=False)
print(summary.to_string(index=False))
print("\nTOP CONCORDANT DOWN")
print(z[z.direction_concordant&(z.meta_effect<0)].sort_values("meta_fdr")[["gene","effect_paired","effect_external","meta_effect","meta_ci_low","meta_ci_high","meta_fdr","I2_percent"]].head(30).to_string(index=False))
