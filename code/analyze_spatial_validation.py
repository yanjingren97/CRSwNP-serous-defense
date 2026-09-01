from pathlib import Path
import re
import sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
xlsx = ROOT / "data/raw/spatial/Validation_cohort_Q3NormalizationFile.xlsx"
module_file = ROOT / "results/locked_40_gene_module.csv"
outdir = ROOT / "results/spatial_validation"
outdir.mkdir(parents=True, exist_ok=True)

mat = pd.read_excel(xlsx, sheet_name="TargetCountMatrix", index_col=0)
module = pd.read_csv(module_file)["gene"].astype(str).tolist()
overlap = [g for g in module if g in mat.index]

# Within-compartment gene standardization prevents highly expressed genes dominating.
rows = []
for compartment in ("PanCK", "CD45"):
    cols = [c for c in mat.columns if c.endswith("_" + compartment)]
    x = np.log2(mat.loc[overlap, cols].astype(float) + 1.0)
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    scores = z.mean(axis=0)
    for sample, score in scores.items():
        m = re.match(r"(CRSwNP|CTRL)_(\d+)_\d+_(PanCK|CD45)$", sample)
        if not m:
            continue
        rows.append({"sample": sample, "group": m.group(1), "patient": f"{m.group(1)}_{m.group(2)}",
                     "compartment": compartment, "score": float(score)})

roi = pd.DataFrame(rows)
patient = roi.groupby(["group", "patient", "compartment"], as_index=False)["score"].mean()

def hedges_g(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n1, n2 = len(a), len(b)
    sp = np.sqrt(((n1-1)*a.var(ddof=1)+(n2-1)*b.var(ddof=1))/(n1+n2-2))
    d = (a.mean()-b.mean())/sp
    return d * (1 - 3/(4*(n1+n2)-9))

summary = []
for compartment in ("PanCK", "CD45"):
    a = patient.query("group == 'CRSwNP' and compartment == @compartment")["score"].to_numpy()
    b = patient.query("group == 'CTRL' and compartment == @compartment")["score"].to_numpy()
    u = mannwhitneyu(a, b, alternative="two-sided", method="exact")
    summary.append({"compartment": compartment, "n_crs": len(a), "n_ctrl": len(b),
                    "mean_crs": a.mean(), "mean_ctrl": b.mean(), "hedges_g_crs_minus_ctrl": hedges_g(a,b),
                    "mannwhitney_u": u.statistic, "p_exact": u.pvalue})

pd.DataFrame({"gene": module, "available": [g in overlap for g in module]}).to_csv(outdir/"module_gene_overlap.csv", index=False)
roi.to_csv(outdir/"roi_module_scores.csv", index=False)
patient.to_csv(outdir/"patient_module_scores.csv", index=False)
pd.DataFrame(summary).to_csv(outdir/"spatial_module_summary.csv", index=False)
print(f"Module overlap: {len(overlap)}/{len(module)}")
print(pd.DataFrame(summary).to_string(index=False))
print("\nPanCK patient scores")
print(patient.query("compartment == 'PanCK'").sort_values(["group","patient"]).to_string(index=False))
