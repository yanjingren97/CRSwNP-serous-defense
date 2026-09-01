from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from secretory_subtype_audit import PANELS, EPI, IMMUNE

BASE = ROOT / "results" / "secretory_subtypes"
DEFINITION_GENES = set(sum(PANELS.values(), []) + EPI + IMMUNE)


def load(ds):
    x = pd.read_csv(BASE/f"{ds}_subtype_counts.tsv.gz", sep="\t", index_col=0)
    m = pd.read_csv(BASE/f"{ds}_subtype_metadata.csv").set_index("column").loc[x.columns]
    lib = x.sum(axis=0).replace(0, np.nan)
    return np.log2(x.div(lib, axis=1)*1e6 + .5), m


def unpaired_vector(x, m, subtype, tissue_a, tissue_b, min_cells):
    z = m[(m.disease == "CRSwNP") & (m.subtype == subtype) &
          m.tissue.isin([tissue_a, tissue_b]) & (m.n_cells >= min_cells)]
    a, b = z[z.tissue == tissue_a], z[z.tissue == tissue_b]
    if len(a) < 2 or len(b) < 2: return None
    keep = ((x[a.index] > 0).sum(axis=1) + (x[b.index] > 0).sum(axis=1)) >= 2
    y = x.loc[keep]
    return pd.DataFrame({"gene": y.index, "effect": y[a.index].mean(axis=1)-y[b.index].mean(axis=1),
                         "n_a": len(a), "n_b": len(b)})


def paired_vector(x, m, subtype, tissue_a, tissue_b, min_cells):
    z = m[(m.disease == "CRSwNP") & (m.subtype == subtype) &
          m.tissue.isin([tissue_a, tissue_b]) & (m.n_cells >= min_cells)].copy()
    z["column_id"] = z.index
    w = z.reset_index(drop=True).pivot(index="donor", columns="tissue", values="column_id").dropna()
    if tissue_a not in w or tissue_b not in w or len(w) < 2: return None
    a, b = w[tissue_a].tolist(), w[tissue_b].tolist()
    keep = ((x[a] > 0).sum(axis=1) + (x[b] > 0).sum(axis=1)) >= 2
    y = x.loc[keep]
    delta = y[a].to_numpy()-y[b].to_numpy()
    return pd.DataFrame({"gene": y.index, "effect": delta.mean(axis=1), "n_a": len(a), "n_b": len(b),
                         "pair_direction_consistency": np.maximum((delta > 0).mean(axis=1), (delta < 0).mean(axis=1))})


def leave_one_out_stability(z):
    # Correlation of the full vector with itself after deleting each paired donor is
    # computed upstream only when raw paired matrices are available; here report
    # robust non-definition-gene concordance between studies.
    strong = z[(z.effect_paired.abs() >= .5) & (z.effect_external.abs() >= .5)]
    return len(strong), ((np.sign(strong.effect_paired) == np.sign(strong.effect_external)).mean() if len(strong) else np.nan)


def main():
    x235, m235 = load("GSE235711")
    x276, m276 = load("GSE276503")
    effects, summary = [], []
    targets = ["Serous_glandular", "Mucous_glandular", "Inflammatory_secretory", "Transitional_secretory"]
    for threshold in [10, 20, 50]:
        for st in targets:
            p = paired_vector(x235, m235, st, "Nasal_polyp", "Ethmoid", threshold)
            u = unpaired_vector(x276, m276, st, "Nasal_polyp", "Inferior_turbinate", threshold)
            if p is None or u is None: continue
            p2 = p[~p.gene.isin(DEFINITION_GENES)]
            u2 = u[~u.gene.isin(DEFINITION_GENES)]
            z = p2[["gene", "effect"]].merge(u2[["gene", "effect"]], on="gene", suffixes=("_paired", "_external"))
            rho, pv = spearmanr(z.effect_paired, z.effect_external)
            strong, rate = leave_one_out_stability(z)
            summary.append({"cell_min": threshold, "subtype": st, "n_paired": int(p.n_a.iloc[0]),
                            "n_external_polyp": int(u.n_a.iloc[0]), "n_external_IT": int(u.n_b.iloc[0]),
                            "genes_non_definition": len(z), "rho": rho, "p": pv,
                            "strong_non_definition": strong, "same_direction_rate": rate})
            z.insert(0, "subtype", st); z.insert(0, "cell_min", threshold)
            effects.append(z)
    summary = pd.DataFrame(summary)
    summary.to_csv(BASE/"subtype_nondefinition_concordance.csv", index=False)
    if effects:
        pd.concat(effects, ignore_index=True).to_csv(BASE/"subtype_nondefinition_gene_effects.csv.gz", index=False, compression="gzip")
    print(summary.to_string(index=False))
    if effects:
        q = pd.concat(effects, ignore_index=True)
        q = q[(q.cell_min == 20) & (q.effect_paired.abs() >= 1) & (q.effect_external.abs() >= 1) &
              (np.sign(q.effect_paired) == np.sign(q.effect_external))]
        q["min_abs"] = q[["effect_paired", "effect_external"]].abs().min(axis=1)
        print("\nTOP NON-DEFINITION GENES")
        print(q.sort_values(["subtype", "min_abs"], ascending=[True, False]).groupby("subtype").head(20).to_string(index=False))


if __name__ == "__main__": main()
