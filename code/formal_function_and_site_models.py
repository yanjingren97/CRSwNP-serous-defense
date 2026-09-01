from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, t as tdist

OUT = ROOT / "results" / "formal_models"
OUT.mkdir(parents=True, exist_ok=True)

def bh(p):
    p=np.asarray(p,float); n=len(p); order=np.argsort(p); q=np.empty(n); last=1.0
    for rank, idx in reversed(list(enumerate(order,1))):
        last=min(last,p[idx]*n/rank); q[idx]=last
    return q

def read_gmt(path):
    sets={}
    with open(path,encoding="utf-8",errors="replace") as h:
        for line in h:
            a=line.rstrip("\n").split("\t")
            if len(a)>=3: sets[a[0]]=set(a[2:])
    return sets

def enrich(query, universe, sets, library, direction):
    q=set(query)&set(universe); u=set(universe); rows=[]
    for term, genes in sets.items():
        g=genes&u; hit=q&g
        if len(g)<10 or len(g)>1000 or len(hit)<3: continue
        p=hypergeom.sf(len(hit)-1,len(u),len(g),len(q))
        rows.append({"library":library,"direction":direction,"term":term,"overlap":len(hit),
                     "query_size":len(q),"term_size":len(g),"p":p,"genes":";".join(sorted(hit))})
    z=pd.DataFrame(rows)
    if len(z): z["fdr"]=bh(z.p); z=z.sort_values(["fdr","p","overlap"])
    return z

effects=pd.read_csv(ROOT/"results/secretory_subtypes/subtype_nondefinition_gene_effects.csv.gz")
e=effects[(effects.cell_min==20)&(effects.subtype=="Serous_glandular")].drop_duplicates("gene")
universe=set(e.gene)
down=e[(e.effect_paired<=-.5)&(e.effect_external<=-.5)].gene
up=e[(e.effect_paired>=.5)&(e.effect_external>=.5)].gene
all_enr=[]
for lib,fn in [("GO_BP_2025","GO_Biological_Process_2025.txt"),("Reactome_2024","Reactome_Pathways_2024.txt")]:
    gs=read_gmt(ROOT/"data/resources/gene_sets"/fn)
    all_enr += [enrich(down,universe,gs,lib,"down"),enrich(up,universe,gs,lib,"up")]
enr=pd.concat(all_enr,ignore_index=True)
enr.to_csv(OUT/"serous_concordant_functional_enrichment.csv",index=False)
pd.DataFrame({"direction":["down","up"],"n_genes":[len(down),len(up)],"background":[len(universe)]*2}).to_csv(OUT/"enrichment_input_summary.csv",index=False)

def contrast_curvature(path,dataset,groups):
    d=pd.read_csv(path); vals={g:d.loc[d.group==g,"score"].to_numpy() for g in groups}
    # (polyp - nonpolyp CRSwNP) - (nonpolyp CRSwNP - healthy): a - 2b + c
    a,b,c=[vals[g] for g in groups]; co=np.array([1.,-2.,1.]); arr=[a,b,c]
    est=a.mean()-2*b.mean()+c.mean()
    var=sum((co[i]**2)*arr[i].var(ddof=1)/len(arr[i]) for i in range(3))
    se=np.sqrt(var); df=var**2/sum(((co[i]**2*arr[i].var(ddof=1)/len(arr[i]))**2)/(len(arr[i])-1) for i in range(3))
    tv=est/se; p=2*tdist.sf(abs(tv),df)
    return {"dataset":dataset,"contrast":"(NP-CRSwNP_nonpolyp)-(CRSwNP_nonpolyp-healthy)",
            "n_np":len(a),"n_crs_nonpolyp":len(b),"n_healthy":len(c),"estimate":est,"se":se,"t":tv,"df":df,"p":p,
            "mean_np":a.mean(),"mean_crs_nonpolyp":b.mean(),"mean_healthy":c.mean()}

models=pd.DataFrame([
    contrast_curvature(ROOT/"results/bulk_crs_validation/GSE36830_scores.csv","GSE36830",["CRSwNP_NP","CRSwNP_UT","Healthy_UT"]),
    contrast_curvature(ROOT/"results/bulk_crs_validation/GSE136825_scores.csv","GSE136825",["CRSwNP_NP","CRSwNP_IT","Healthy_IT"])
])
models["fdr_two_tests"]=bh(models.p)
models.to_csv(OUT/"anatomy_disease_contrast_of_contrasts.csv",index=False)
print("ENRICHMENT INPUT",len(down),"down",len(up),"up; background",len(universe))
print(enr[enr.fdr<.05].groupby(["library","direction"]).head(12)[["library","direction","term","overlap","term_size","p","fdr"]].to_string(index=False))
print("\nCONTRAST OF CONTRASTS")
print(models.to_string(index=False))
