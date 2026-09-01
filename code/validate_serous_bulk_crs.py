from pathlib import Path
from io import StringIO
import csv, gzip, sys
ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

RAW = ROOT / "data" / "raw" / "bulk_crs"
OUT = ROOT / "results" / "bulk_crs_validation"
OUT.mkdir(parents=True, exist_ok=True)


def hedges(a, b):
    a, b = np.asarray(a), np.asarray(b); n1, n0 = len(a), len(b)
    sp = np.sqrt(((n1-1)*a.var(ddof=1)+(n0-1)*b.var(ddof=1))/(n1+n0-2))
    if sp == 0: return np.nan
    return ((a.mean()-b.mean())/sp)*(1-3/(4*(n1+n0)-9))


def read_geo(path):
    meta, lines, inside = {}, [], False
    with gzip.open(path, "rt", errors="replace") as h:
        for line in h:
            if line.startswith("!Sample_"):
                parts = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
                meta.setdefault(parts[0], []).append(parts[1:])
            if line.startswith("!series_matrix_table_begin"): inside = True; continue
            if line.startswith("!series_matrix_table_end"): inside = False; continue
            if inside: lines.append(line)
    mat = pd.read_csv(StringIO("".join(lines)), sep="\t", index_col=0)
    return mat, meta


def annotation(path):
    with gzip.open(path, "rt", errors="replace") as h:
        lines = [line for line in h if not line.startswith(("#", "^", "!"))]
    a = pd.read_csv(StringIO("".join(lines)), sep="\t", dtype=str)
    symbol_col = next(c for c in a.columns if c.lower().replace("_", " ") in ["gene symbol", "gene assignment"] or c.lower()=="gene symbol")
    return a.set_index("ID")[symbol_col].fillna("").str.split(" /// ").str[0].str.split("//").str[0].str.strip()


def collapse(expr, symbols):
    z = expr.copy(); z["symbol"] = symbols.reindex(z.index).fillna("").values
    z = z[z.symbol.ne("")].groupby("symbol").mean()
    return z


def score(expr, module):
    g = [x for x in module if x in expr.index]
    x = expr.loc[g]
    sd = x.std(axis=1).replace(0, np.nan)
    z = x.sub(x.mean(axis=1), axis=0).div(sd, axis=0)
    return z.mean(axis=0), g


def contrasts(dataset, scores, groups, pairs):
    rows=[]
    for case, ctrl in pairs:
        a=scores[groups==case]; b=scores[groups==ctrl]
        if len(a)<2 or len(b)<2: continue
        rows.append({"dataset":dataset,"case":case,"control":ctrl,"n_case":len(a),"n_control":len(b),
                     "median_case":np.median(a),"median_control":np.median(b),
                     "median_difference":np.median(a)-np.median(b),"hedges_g":hedges(a,b),
                     "mannwhitney_p":mannwhitneyu(a,b,alternative="two-sided").pvalue})
    return rows


def module_genes():
    e=pd.read_csv(ROOT/"results"/"secretory_subtypes"/"subtype_nondefinition_gene_effects.csv.gz")
    z=e[(e.cell_min==20)&(e.subtype=="Serous_glandular")&(e.effect_paired<=-1)&(e.effect_external<=-1)].copy()
    z["min_abs"]=z[["effect_paired","effect_external"]].abs().min(axis=1)
    banned={"AQP1","ACKR1","EMCN","VWF","PECAM1","KDR","COL1A1","COL1A2","COL3A1","DCN","LUM","HBB","HBA1","HBA2"}
    z=z[~z.gene.isin(banned)&~z.gene.str.startswith(("IG","HLA-","RPL","RPS"))]
    z=z[~z.gene.str.contains(r"\.",regex=True)]
    genes=z.sort_values("min_abs",ascending=False).gene.drop_duplicates().head(40).tolist()
    pd.DataFrame({"gene":genes}).to_csv(ROOT/"results/locked_40_gene_module.csv",index=False)
    return genes


def gse36830(module):
    x,m=read_geo(RAW/"GSE36830_series_matrix.txt.gz")
    sy=annotation(RAW/"GPL570.annot.gz")
    x=collapse(x,sy)
    sources=np.array(m["!Sample_source_name_ch1"][0])
    titles=np.array(m["!Sample_title"][0])
    groups=[]
    for s,t in zip(sources,titles):
        text=(s+" "+t).lower()
        if "control" in text and "uncinate" in text: groups.append("Healthy_UT")
        elif "without nasal polyps" in text or "crssnp" in text: groups.append("CRSsNP_UT")
        elif "nasal polyp tissue" in text or ("polyp" in text and "uncinate" not in text): groups.append("CRSwNP_NP")
        elif "with nasal polyps" in text or "crswnp" in text: groups.append("CRSwNP_UT")
        else: groups.append(t)
    s,g=score(x,module); groups=np.array(groups)
    pd.DataFrame({"sample":x.columns,"title":titles,"source":sources,"group":groups,"score":s.values}).to_csv(OUT/"GSE36830_scores.csv",index=False)
    return contrasts("GSE36830",s,groups,[("CRSwNP_NP","Healthy_UT"),("CRSwNP_UT","Healthy_UT"),("CRSwNP_NP","CRSwNP_UT"),("CRSsNP_UT","Healthy_UT")]),g


def gse136825(module):
    counts=pd.read_csv(RAW/"GSE136825_genecounts_20190903.txt.gz",sep="\t",index_col=0)
    # Geneid is already the index; five remaining featureCounts annotation columns precede samples.
    sample_cols=counts.columns[5:]
    counts=counts[sample_cols]
    # Map Ensembl IDs using the 10x feature table already downloaded for discovery data.
    manifest=pd.read_csv(ROOT/"results"/"preflight"/"sample_manifest.csv")
    feature_path=manifest.loc[manifest.dataset=="GSE235711","features"].iloc[0]
    feat=pd.read_csv(feature_path,sep="\t",header=None)
    mp=pd.Series(feat.iloc[:,1].values,index=feat.iloc[:,0].str.split(".").str[0]).drop_duplicates()
    counts.index=counts.index.str.split(".").str[0]
    counts["symbol"]=mp.reindex(counts.index).fillna("").values
    counts=counts[counts.symbol.ne("")].groupby("symbol").sum()
    lib=counts.sum(axis=0).replace(0,np.nan)
    expr=np.log2(counts.div(lib,axis=1)*1e6+.5)
    _,m=read_geo(RAW/"GSE136825_series_matrix.txt.gz")
    sources=np.array(m["!Sample_source_name_ch1"][0])
    titles=np.array(m["!Sample_title"][0])
    if len(sources)!=expr.shape[1]: raise ValueError(f"metadata/count mismatch {len(sources)} vs {expr.shape[1]}")
    groups=np.array(["CRSwNP_NP" if s=="Nasal Polyp Tissue" else "CRSwNP_IT" if s=="Nasal Polyp Inferior Turbinate" else "Healthy_IT" for s in sources])
    s,g=score(expr,module)
    pd.DataFrame({"sample":expr.columns,"title":titles,"source":sources,"group":groups,"score":s.values}).to_csv(OUT/"GSE136825_scores.csv",index=False)
    return contrasts("GSE136825",s,groups,[("CRSwNP_NP","Healthy_IT"),("CRSwNP_IT","Healthy_IT"),("CRSwNP_NP","CRSwNP_IT")]),g


def main():
    module=module_genes(); rows=[]
    a,g1=gse36830(module); rows+=a
    b,g2=gse136825(module); rows+=b
    out=pd.DataFrame(rows); out.to_csv(OUT/"serous_module_validation.csv",index=False)
    print("MODULE",module); print("genes GSE36830",len(g1),"genes GSE136825",len(g2)); print(out.to_string(index=False))


if __name__=="__main__": main()
