from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

BASE = ROOT / "results" / "celltype_preflight"


def hedges(case, ctrl):
    n1, n0 = len(case), len(ctrl)
    if n1 < 2 or n0 < 2:
        return np.nan
    sp = np.sqrt(((n1 - 1) * np.var(case, ddof=1) + (n0 - 1) * np.var(ctrl, ddof=1)) / (n1 + n0 - 2))
    if sp == 0 or not np.isfinite(sp): return np.nan
    return ((np.mean(case) - np.mean(ctrl)) / sp) * (1 - 3 / (4 * (n1 + n0) - 9))


def compare(df, value, ds, tissue, case, ctrl, feature_col, features):
    out = []
    sub = df[(df.dataset == ds) & (df.tissue == tissue)]
    for feat in features:
        q = sub[sub[feature_col] == feat]
        a = q[q.disease == case][value].dropna().to_numpy()
        b = q[q.disease == ctrl][value].dropna().to_numpy()
        if len(a) and len(b):
            out.append({"dataset": ds, "tissue": tissue, "case": case, "control": ctrl,
                        "feature": feat, "n_case": len(a), "n_control": len(b),
                        "case_median": np.median(a), "control_median": np.median(b),
                        "median_difference": np.median(a)-np.median(b), "hedges_g": hedges(a,b),
                        "mannwhitney_p": mannwhitneyu(a,b,alternative="two-sided").pvalue})
    return out


def main():
    comp = pd.read_csv(BASE / "donor_celltype_composition.csv")
    epi = comp[comp.level == "epithelial"]
    features = ["Basal", "Ciliated", "Secretory", "Goblet", "Cycling", "Ionocyte", "Tuft"]
    contrasts = [("GSE235711","Ethmoid","CRSsNP","Healthy"),
                 ("GSE235711","Ethmoid","CRSwNP","Healthy"),
                 ("GSE276503","Inferior_turbinate","CRSwNP","Healthy"),
                 ("GSE261706","Inferior_turbinate","AR","Healthy")]
    out=[]
    for c in contrasts:
        out += compare(epi,"proportion",*c,"celltype",features)
    pd.DataFrame(out).to_csv(BASE / "patient_epithelial_composition_effects.csv",index=False)

    # Aggregate technical libraries, calculate compartment-specific marker CPM.
    p = pd.read_csv(BASE / "sample_compartment_signatures.csv")
    keys=["dataset","donor","disease","tissue","compartment","signature"]
    p=p.groupby(keys,as_index=False).agg(n_cells=("n_cells","sum"),marker_counts=("marker_counts","sum"),
                                         all_counts=("all_counts","sum"),genes_present=("genes_present","max"))
    p["score"] = np.log1p(1e6*p.marker_counts/p.all_counts.clip(lower=1)/p.genes_present.clip(lower=1))
    sigs=p.signature.unique().tolist()
    sout=[]
    for c in contrasts:
        for compartment in ["Epithelial","Fibroblast"]:
            q=p[(p.compartment==compartment)&(p.n_cells>=50)]
            z=compare(q,"score",*c,"signature",sigs)
            for item in z: item["compartment"]=compartment
            sout+=z
    pd.DataFrame(sout).to_csv(BASE / "patient_compartment_signature_effects.csv",index=False)

    # Within-patient polyp-vs-ethmoid effects in the four explicitly paired donors.
    paired=[]
    q=epi[(epi.dataset=="GSE235711")&(epi.disease=="CRSwNP")]
    for feat in features:
        w=q[q.celltype==feat].pivot(index="donor",columns="tissue",values="proportion").dropna()
        if len(w):
            diff=w.Nasal_polyp-w.Ethmoid
            paired.append({"feature":feat,"n_pairs":len(w),"median_paired_difference":diff.median(),
                           "positive_pairs":int((diff>0).sum()),
                           "wilcoxon_p":wilcoxon(diff).pvalue if len(w)>1 and np.any(diff!=0) else np.nan})
    pd.DataFrame(paired).to_csv(BASE / "paired_polyp_ethmoid_composition.csv",index=False)
    print(pd.DataFrame(out).query("feature in ['Basal','Ciliated','Goblet','Secretory']").to_string(index=False))
    print("\nPAIRED\n",pd.DataFrame(paired).to_string(index=False))


if __name__ == "__main__": main()
