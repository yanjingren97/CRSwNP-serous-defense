from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.special import logit
from scipy.stats import mannwhitneyu, wilcoxon

OUT=ROOT/"results/frozen_v1"; OUT.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(20260831)

def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); last=1.
    for rank,i in reversed(list(enumerate(o,1))): last=min(last,p[i]*n/rank); q[i]=last
    return q

def hedges(a,b):
    a,b=np.asarray(a),np.asarray(b); n1,n0=len(a),len(b)
    sp=np.sqrt(((n1-1)*a.var(ddof=1)+(n0-1)*b.var(ddof=1))/(n1+n0-2))
    return ((a.mean()-b.mean())/sp)*(1-3/(4*(n1+n0)-9)) if sp>0 else np.nan

def ci_unpaired(a,b,n=10000):
    vals=np.array([np.median(rng.choice(a,len(a),True))-np.median(rng.choice(b,len(b),True)) for _ in range(n)])
    return np.quantile(vals,[.025,.975])

def ci_paired(d,n=10000):
    vals=np.array([np.median(rng.choice(d,len(d),True)) for _ in range(n)])
    return np.quantile(vals,[.025,.975])

# Freeze sample inclusion/exclusion ledger.
manifest=pd.read_csv(ROOT/"results/preflight/sample_manifest.csv")
errors=pd.read_csv(ROOT/"results/secretory_subtypes/errors.csv")
ledger=manifest[["dataset","sample","donor","disease","tissue","library"]].copy()
ledger=ledger.merge(errors,on=["dataset","sample"],how="left")
ledger["included_readable"]=ledger.error.isna()
ledger["exclusion_reason"]=ledger.error.fillna("")
ledger.to_csv(OUT/"sample_inclusion_ledger.csv",index=False)

d=pd.read_csv(ROOT/"results/secretory_subtypes/donor_subtype_composition.csv")
d["classification_coverage"]=d.n_classified_secretory/d.n_secretory_candidates.clip(lower=1)
d["logit_candidate_fraction"]=logit((d.n_cells+.5)/(d.n_secretory_candidates+1))
d.to_csv(OUT/"patient_secretory_composition.csv",index=False)

# Freeze all Serous pseudobulk columns and primary eligibility at 10/20/50 cells.
meta=[]
for ds in ["GSE235711","GSE276503","GSE261706"]:
    z=pd.read_csv(ROOT/f"results/secretory_subtypes/{ds}_subtype_metadata.csv")
    z=z[z.subtype=="Serous_glandular"].copy()
    for k in [10,20,50]: z[f"eligible_n{k}"]=z.n_cells>=k
    meta.append(z)
pd.concat(meta,ignore_index=True).to_csv(OUT/"serous_pseudobulk_metadata.csv",index=False)

rows=[]
contrasts=[
 ("unpaired","GSE235711","Ethmoid","CRSsNP","Healthy","CRSsNP_vs_Healthy_Ethmoid"),
 ("unpaired","GSE235711","Ethmoid","CRSwNP","Healthy","CRSwNP_vs_Healthy_Ethmoid"),
 ("unpaired","GSE276503","Inferior_turbinate","CRSwNP","Healthy","CRSwNP_vs_Healthy_IT"),
 ("unpaired","GSE261706","Inferior_turbinate","AR","Healthy","AR_vs_Healthy_IT"),
 ("unpaired_tissue","GSE276503","CRSwNP","Nasal_polyp","Inferior_turbinate","NP_vs_CRSwNP_IT"),
 ("paired","GSE235711","CRSwNP","Nasal_polyp","Ethmoid","paired_NP_vs_Ethmoid")]

for typ,ds,a1,a2,a3,label in contrasts:
  for st in sorted(d.subtype.unique()):
    if typ=="unpaired":
      z=d[(d.dataset==ds)&(d.tissue==a1)&d.disease.isin([a2,a3])&(d.subtype==st)]
      a=z[z.disease==a2]; b=z[z.disease==a3]
    elif typ=="unpaired_tissue":
      z=d[(d.dataset==ds)&(d.disease==a1)&d.tissue.isin([a2,a3])&(d.subtype==st)]
      a=z[z.tissue==a2]; b=z[z.tissue==a3]
    else:
      z=d[(d.dataset==ds)&(d.disease==a1)&d.tissue.isin([a2,a3])&(d.subtype==st)]
      w=z.pivot(index="donor",columns="tissue",values=["proportion_of_candidates","logit_candidate_fraction"]).dropna()
      raw=(w[("proportion_of_candidates",a2)]-w[("proportion_of_candidates",a3)]).to_numpy()
      ld=(w[("logit_candidate_fraction",a2)]-w[("logit_candidate_fraction",a3)]).to_numpy()
      if len(raw)<2: continue
      try: p=wilcoxon(ld,method="auto").pvalue
      except ValueError: p=1.
      lo,hi=ci_paired(raw)
      rows.append({"dataset":ds,"contrast":label,"subtype":st,"design":"paired","n_case":len(raw),"n_control":len(raw),
                   "median_candidate_fraction_case":w[("proportion_of_candidates",a2)].median(),
                   "median_candidate_fraction_control":w[("proportion_of_candidates",a3)].median(),
                   "median_difference":np.median(raw),"ci95_low":lo,"ci95_high":hi,"effect_logit_mean_delta":ld.mean(),
                   "hedges_g_logit":np.nan,"p":p,"positive_pairs":int((raw>0).sum()),"negative_pairs":int((raw<0).sum())})
      continue
    if len(a)<2 or len(b)<2: continue
    av=a.proportion_of_candidates.to_numpy(); bv=b.proportion_of_candidates.to_numpy()
    lo,hi=ci_unpaired(av,bv)
    p=mannwhitneyu(a.logit_candidate_fraction,b.logit_candidate_fraction,alternative="two-sided",method="auto").pvalue
    rows.append({"dataset":ds,"contrast":label,"subtype":st,"design":"unpaired","n_case":len(a),"n_control":len(b),
                 "median_candidate_fraction_case":np.median(av),"median_candidate_fraction_control":np.median(bv),
                 "median_difference":np.median(av)-np.median(bv),"ci95_low":lo,"ci95_high":hi,
                 "effect_logit_mean_delta":a.logit_candidate_fraction.mean()-b.logit_candidate_fraction.mean(),
                 "hedges_g_logit":hedges(a.logit_candidate_fraction,b.logit_candidate_fraction),"p":p,
                 "positive_pairs":np.nan,"negative_pairs":np.nan})

effects=pd.DataFrame(rows)
effects["fdr_within_contrast"]=effects.groupby("contrast").p.transform(lambda x:bh(x.to_numpy()))
effects["endpoint_role"]=np.where(effects.subtype=="Serous_glandular","primary_serous","secondary_subtype")
effects.to_csv(OUT/"formal_composition_effects.csv",index=False)

overview=ledger.groupby(["dataset","disease","tissue"],as_index=False).agg(
    libraries_total=("sample","size"),libraries_included=("included_readable","sum"),donors=("donor","nunique"))
pc=d.groupby(["dataset","disease","tissue"],as_index=False).agg(secretory_candidates=("n_secretory_candidates","sum"),
    classified_secretory=("n_classified_secretory","sum"))
overview.merge(pc,on=["dataset","disease","tissue"],how="left").to_csv(OUT/"cohort_overview.csv",index=False)

print("COHORT OVERVIEW")
print(overview.to_string(index=False))
print("\nPRIMARY SEROUS COMPOSITION")
print(effects[effects.subtype=="Serous_glandular"].to_string(index=False))
