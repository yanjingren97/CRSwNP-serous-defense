from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

INP = ROOT / "results" / "secretory_subtypes" / "donor_subtype_composition.csv"
OUT = ROOT / "results" / "secretory_subtypes"


def independent(d, dataset, tissue, case, control):
    z = d[(d.dataset == dataset) & (d.tissue == tissue) & d.disease.isin([case, control])]
    ans = []
    for st in sorted(z.subtype.unique()):
        a = z[(z.disease == case) & (z.subtype == st)].proportion.to_numpy()
        b = z[(z.disease == control) & (z.subtype == st)].proportion.to_numpy()
        if len(a) < 2 or len(b) < 2: continue
        _, p = mannwhitneyu(a, b, alternative="two-sided", method="auto")
        ans.append({"dataset": dataset, "contrast": f"{case}_vs_{control}_{tissue}", "subtype": st,
                    "n_case": len(a), "n_control": len(b), "median_case": np.median(a),
                    "median_control": np.median(b), "median_difference": np.median(a)-np.median(b), "p": p})
    return ans


def paired(d):
    z = d[(d.dataset == "GSE235711") & (d.disease == "CRSwNP") &
          d.tissue.isin(["Nasal_polyp", "Ethmoid"])]
    ans = []
    for st in sorted(z.subtype.unique()):
        w = z[z.subtype == st].pivot(index="donor", columns="tissue", values="proportion").dropna()
        delta = w.Nasal_polyp - w.Ethmoid
        if len(delta) < 2: continue
        try: p = wilcoxon(delta).pvalue
        except ValueError: p = 1.0
        ans.append({"dataset": "GSE235711", "contrast": "paired_polyp_vs_ethmoid", "subtype": st,
                    "n_case": len(delta), "n_control": len(delta), "median_case": np.median(w.Nasal_polyp),
                    "median_control": np.median(w.Ethmoid), "median_difference": np.median(delta), "p": p,
                    "positive_pairs": int((delta > 0).sum()), "negative_pairs": int((delta < 0).sum())})
    return ans


def tissue_independent(d, dataset, disease, tissue_a, tissue_b):
    z = d[(d.dataset == dataset) & (d.disease == disease) & d.tissue.isin([tissue_a, tissue_b])]
    ans = []
    for st in sorted(z.subtype.unique()):
        a = z[(z.tissue == tissue_a) & (z.subtype == st)].proportion.to_numpy()
        b = z[(z.tissue == tissue_b) & (z.subtype == st)].proportion.to_numpy()
        if len(a) < 2 or len(b) < 2: continue
        _, p = mannwhitneyu(a, b, alternative="two-sided", method="auto")
        ans.append({"dataset": dataset, "contrast": f"{tissue_a}_vs_{tissue_b}", "subtype": st,
                    "n_case": len(a), "n_control": len(b), "median_case": np.median(a),
                    "median_control": np.median(b), "median_difference": np.median(a)-np.median(b), "p": p})
    return ans


def main():
    d = pd.read_csv(INP)
    out = []
    out += independent(d, "GSE235711", "Ethmoid", "CRSsNP", "Healthy")
    out += independent(d, "GSE235711", "Ethmoid", "CRSwNP", "Healthy")
    out += independent(d, "GSE276503", "Inferior_turbinate", "CRSwNP", "Healthy")
    out += independent(d, "GSE261706", "Inferior_turbinate", "AR", "Healthy")
    out += paired(d)
    out += tissue_independent(d, "GSE276503", "CRSwNP", "Nasal_polyp", "Inferior_turbinate")
    tab = pd.DataFrame(out)
    tab.to_csv(OUT / "patient_subtype_composition_effects.csv", index=False)
    print(tab.to_string(index=False))


if __name__ == "__main__": main()
