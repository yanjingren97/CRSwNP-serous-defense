from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from validate_serous_bulk_crs import read_geo, annotation, collapse, score, hedges, RAW

OUT=ROOT/"results"/"bulk_crs_validation"; OUT.mkdir(parents=True,exist_ok=True)
ECM=["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","DCN","LUM","FN1","SPARC","TGFBI","POSTN"]


def ordered_module():
    e=pd.read_csv(ROOT/"results"/"secretory_subtypes"/"subtype_nondefinition_gene_effects.csv.gz")
    z=e[(e.cell_min==20)&(e.subtype=="Serous_glandular")&(e.effect_paired<=-1)&(e.effect_external<=-1)].copy()
    z["min_abs"]=z[["effect_paired","effect_external"]].abs().min(axis=1)
    banned={"AQP1","ACKR1","EMCN","VWF","PECAM1","KDR","COL1A1","COL1A2","COL3A1","DCN","LUM","HBB","HBA1","HBA2"}
    z=z[~z.gene.isin(banned)&~z.gene.str.startswith(("IG","HLA-","RPL","RPS"))&~z.gene.str.contains(r"\.",regex=True)]
    return z.sort_values("min_abs",ascending=False).gene.drop_duplicates().head(40).tolist()


def load368():
    x,m=read_geo(RAW/"GSE36830_series_matrix.txt.gz"); x=collapse(x,annotation(RAW/"GPL570.annot.gz"))
    src=np.array(m["!Sample_source_name_ch1"][0]); titles=np.array(m["!Sample_title"][0]); groups=[]
    for s,t in zip(src,titles):
        q=(s+" "+t).lower()
        if "control" in q and "uncinate" in q: groups.append("Healthy_UT")
        elif "without nasal polyps" in q or "crssnp" in q: groups.append("CRSsNP_UT")
        elif "nasal polyp tissue" in q or ("polyp" in q and "uncinate" not in q): groups.append("CRSwNP_NP")
        elif "with nasal polyps" in q or "crswnp" in q: groups.append("CRSwNP_UT")
        else: groups.append(t)
    return x,np.array(groups)


def load136():
    c=pd.read_csv(RAW/"GSE136825_genecounts_20190903.txt.gz",sep="\t",index_col=0); c=c[c.columns[5:]]
    man=pd.read_csv(ROOT/"results"/"preflight"/"sample_manifest.csv"); fp=man.loc[man.dataset=="GSE235711","features"].iloc[0]
    f=pd.read_csv(fp,sep="\t",header=None); mp=pd.Series(f.iloc[:,1].values,index=f.iloc[:,0].str.split(".").str[0]).drop_duplicates()
    c.index=c.index.str.split(".").str[0]; sy=mp.reindex(c.index).fillna("").values
    c=pd.concat([c.reset_index(drop=True),pd.Series(sy,name="symbol")],axis=1); c=c[c.symbol.ne("")].groupby("symbol").sum()
    x=np.log2(c.div(c.sum(axis=0),axis=1)*1e6+.5)
    _,m=read_geo(RAW/"GSE136825_series_matrix.txt.gz"); src=np.array(m["!Sample_source_name_ch1"][0])
    groups=np.array(["CRSwNP_NP" if s=="Nasal Polyp Tissue" else "CRSwNP_IT" if s=="Nasal Polyp Inferior Turbinate" else "Healthy_IT" for s in src])
    return x,groups


def compare(dataset,x,groups,module,case,ctrl,label):
    s,g=score(x,module); a=s[groups==case]; b=s[groups==ctrl]
    return {"dataset":dataset,"module":label,"genes_detected":len(g),"case":case,"control":ctrl,
            "n_case":len(a),"n_control":len(b),"hedges_g":hedges(a,b),
            "median_difference":np.median(a)-np.median(b),"p":mannwhitneyu(a,b).pvalue}


def main():
    module=ordered_module(); x368,g368=load368(); x136,g136=load136(); rows=[]
    variants={"top10":module[:10],"top20":module[:20],"top40":module,"drop_top5":module[5:]}
    for name,genes in variants.items():
        rows.append(compare("GSE36830",x368,g368,genes,"CRSwNP_NP","CRSwNP_UT",name))
        rows.append(compare("GSE136825",x136,g136,genes,"CRSwNP_NP","CRSwNP_IT",name))
    # Leave-one-gene-out distributions.
    for gene in module:
        genes=[g for g in module if g!=gene]
        r=compare("GSE36830",x368,g368,genes,"CRSwNP_NP","CRSwNP_UT",f"LOO:{gene}"); rows.append(r)
        r=compare("GSE136825",x136,g136,genes,"CRSwNP_NP","CRSwNP_IT",f"LOO:{gene}"); rows.append(r)
    out=pd.DataFrame(rows); out.to_csv(OUT/"serous_module_sensitivity.csv",index=False)

    coupling=[]
    for ds,x,groups in [("GSE36830",x368,g368),("GSE136825",x136,g136)]:
        ser,_=score(x,module); ecm,eg=score(x,ECM)
        for grp in np.unique(groups):
            mask=groups==grp
            if mask.sum()<5: continue
            rho,p=spearmanr(ser[mask],ecm[mask])
            coupling.append({"dataset":ds,"group":grp,"n":int(mask.sum()),"ecm_genes":len(eg),"rho":rho,"p":p})
        # Remove group means to estimate within-group association over all samples.
        sr=pd.Series(ser.values).groupby(groups).transform(lambda v:v-v.mean())
        er=pd.Series(ecm.values).groupby(groups).transform(lambda v:v-v.mean())
        rho,p=spearmanr(sr,er); coupling.append({"dataset":ds,"group":"within_group_residual","n":len(groups),"ecm_genes":len(eg),"rho":rho,"p":p})
    coupling=pd.DataFrame(coupling); coupling.to_csv(OUT/"serous_ecm_coupling.csv",index=False)
    print("VARIANTS"); print(out[~out.module.str.startswith("LOO:")].to_string(index=False))
    print("\nLOO RANGE"); print(out[out.module.str.startswith("LOO:")].groupby("dataset").agg(g_min=("hedges_g","min"),g_max=("hedges_g","max"),p_min=("p","min"),p_max=("p","max")).to_string())
    print("\nCOUPLING"); print(coupling.to_string(index=False))


if __name__=="__main__":main()
