from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

src=ROOT/"data/raw/proteomics/PXD013330_MQ_results/proteinGroups.txt"
outdir=ROOT/"results/protein_validation"; outdir.mkdir(parents=True,exist_ok=True)
x=pd.read_csv(src,sep="\t",low_memory=False)
x=x[(x["Reverse"]!="+")&(x["Potential contaminant"]!="+")&x["Gene names"].notna()].copy()
x["gene"]=x["Gene names"].str.split(";").str[0]
cols={g:[f"LFQ intensity {g}_{i}" for i in range(1,4)] for g in ["CON","CRS","CRSwNP"]}
for c in sum(cols.values(),[]): x[c]=pd.to_numeric(x[c],errors="coerce").replace(0,np.nan)
x=x.sort_values("Score",ascending=False).drop_duplicates("gene").set_index("gene")
module=pd.read_csv(ROOT/"results/locked_40_gene_module.csv")["gene"].tolist()
det=[g for g in module if g in x.index and x.loc[g,sum(cols.values(),[])].notna().sum()>=2]

rows=[]
for g in det:
    r={"gene":g}
    for grp,cs in cols.items():
        v=np.log2(x.loc[g,cs].dropna().astype(float))
        r[f"detected_{grp}"]=len(v); r[f"median_log2_{grp}"]=v.median() if len(v) else np.nan
    r["log2diff_CRSwNP_minus_CON"]=r["median_log2_CRSwNP"]-r["median_log2_CON"]
    rows.append(r)
genes=pd.DataFrame(rows).sort_values("log2diff_CRSwNP_minus_CON")

# Sample-level score uses only proteins quantified in every one of the nine samples.
complete=[g for g in module if g in x.index and x.loc[g,sum(cols.values(),[])].notna().all()]
mat=np.log2(x.loc[complete,sum(cols.values(),[])].astype(float))
z=mat.sub(mat.mean(axis=1),axis=0).div(mat.std(axis=1,ddof=0),axis=0)
scores=z.mean(axis=0)
sample=[]
for grp,cs in cols.items():
    for c in cs: sample.append({"sample":c.replace("LFQ intensity ",""),"group":grp,"score":scores[c]})
sample=pd.DataFrame(sample)
a=sample.query("group=='CRSwNP'").score; b=sample.query("group=='CON'").score
summary=pd.DataFrame([{"dataset":"PXD013330_DDA","module_genes_detected_any":len(det),
                       "module_genes_complete_9_samples":len(complete),"n_CRSwNP":len(a),"n_CON":len(b),
                       "replicate_unit":"technical replicates of pooled biological samples",
                       "median_CRSwNP":a.median(),"median_CON":b.median(),
                       "median_difference":a.median()-b.median(),
                       "exploratory_mannwhitney_p_not_biological":mannwhitneyu(a,b,method="exact").pvalue,
                       "genes_complete":";".join(complete)}])
genes.to_csv(outdir/"PXD013330_module_proteins.csv",index=False)
sample.to_csv(outdir/"PXD013330_sample_module_scores.csv",index=False)
summary.to_csv(outdir/"PXD013330_module_summary.csv",index=False)
print(summary.to_string(index=False)); print("\nPROTEINS"); print(genes.to_string(index=False))
