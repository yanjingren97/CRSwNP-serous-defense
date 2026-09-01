from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from secretory_subtype_pseudobulk_validation import load, unpaired_vector, DEFINITION_GENES

BASE = ROOT / "results" / "secretory_subtypes"
OUT = ROOT / "results" / "robustness"
OUT.mkdir(parents=True, exist_ok=True)

x235, m235 = load("GSE235711")
x276, m276 = load("GSE276503")
st = "Serous_glandular"
threshold = 20

z = m235[(m235.disease == "CRSwNP") & (m235.subtype == st) &
         m235.tissue.isin(["Nasal_polyp", "Ethmoid"]) & (m235.n_cells >= threshold)].copy()
z["column_id"] = z.index
w = z.reset_index(drop=True).pivot(index="donor", columns="tissue", values="column_id").dropna()
external = unpaired_vector(x276, m276, st, "Nasal_polyp", "Inferior_turbinate", threshold)
external = external[~external.gene.isin(DEFINITION_GENES)][["gene", "effect"]].rename(columns={"effect":"external"})
locked = pd.read_csv(ROOT/"results/locked_40_gene_module.csv")["gene"].tolist()

rows = []
gene_rows = []
for omitted in [None] + w.index.tolist():
    ww = w if omitted is None else w.drop(index=omitted)
    a, b = ww["Nasal_polyp"].tolist(), ww["Ethmoid"].tolist()
    y = x235[~x235.index.isin(DEFINITION_GENES)]
    effect = (y[a].to_numpy() - y[b].to_numpy()).mean(axis=1)
    local = pd.DataFrame({"gene": y.index, "paired": effect})
    merged = local.merge(external, on="gene")
    rho, p = spearmanr(merged.paired, merged.external)
    strong = merged[(merged.paired.abs() >= .5) & (merged.external.abs() >= .5)]
    mod = merged[merged.gene.isin(locked)].copy()
    rows.append({"omitted_donor": "none_full" if omitted is None else omitted,
                 "n_pairs": len(ww), "genes": len(merged), "rho_external": rho,
                 "rho_p": p, "strong_genes": len(strong),
                 "strong_same_direction": (np.sign(strong.paired)==np.sign(strong.external)).mean(),
                 "locked_available": len(mod), "locked_negative_paired": int((mod.paired < 0).sum()),
                 "locked_median_paired_effect": mod.paired.median()})
    mod.insert(0, "omitted_donor", "none_full" if omitted is None else omitted)
    gene_rows.append(mod)

pd.DataFrame(rows).to_csv(OUT/"serous_patient_leave_one_out_summary.csv", index=False)
pd.concat(gene_rows, ignore_index=True).to_csv(OUT/"serous_patient_leave_one_out_genes.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
