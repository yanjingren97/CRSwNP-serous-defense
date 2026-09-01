from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_rel, ttest_ind

BASE = ROOT / "results" / "pseudobulk"
OUT = ROOT / "results" / "anatomy_aware"
OUT.mkdir(parents=True, exist_ok=True)
CELL_MIN = 50


def load(ds):
    x = pd.read_csv(BASE / f"{ds}_pseudobulk_counts.tsv.gz", sep="\t", index_col=0)
    m = pd.read_csv(BASE / f"{ds}_pseudobulk_metadata.csv").set_index("column").loc[x.columns]
    lib = x.sum(axis=0).replace(0, np.nan)
    return np.log2(x.div(lib, axis=1) * 1e6 + 0.5), m


def vector(logcpm, meta, disease, tissue_a, tissue_b, celltype, paired=False, cell_min=CELL_MIN):
    m = meta[(meta.disease == disease) & (meta.celltype == celltype) &
             (meta.tissue.isin([tissue_a, tissue_b])) & (meta.n_cells >= cell_min)]
    if paired:
        mm = m.copy()
        mm["column_id"] = mm.index
        wide = mm.reset_index(drop=True).pivot(index="donor", columns="tissue", values="column_id").dropna()
        if tissue_a not in wide or tissue_b not in wide or len(wide) < 2:
            return None
        a = logcpm[wide[tissue_a].tolist()].to_numpy()
        b = logcpm[wide[tissue_b].tolist()].to_numpy()
        stat, p = ttest_rel(a, b, axis=1, nan_policy="omit")
        delta = a - b
        return pd.DataFrame({"gene": logcpm.index, "effect": delta.mean(axis=1),
                             "p": p, "n_a": len(wide), "n_b": len(wide),
                             "sign_consistency": np.maximum((delta > 0).mean(axis=1),
                                                             (delta < 0).mean(axis=1))})
    ma, mb = m[m.tissue == tissue_a], m[m.tissue == tissue_b]
    if len(ma) < 2 or len(mb) < 2:
        return None
    a, b = logcpm[ma.index].to_numpy(), logcpm[mb.index].to_numpy()
    stat, p = ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    return pd.DataFrame({"gene": logcpm.index, "effect": a.mean(axis=1)-b.mean(axis=1),
                         "p": p, "n_a": len(ma), "n_b": len(mb),
                         "sign_consistency": np.nan})


def main():
    x235, m235 = load("GSE235711")
    x276, m276 = load("GSE276503")
    disease = pd.read_csv(BASE / "all_pseudobulk_gene_effects.csv.gz")
    rows, concord = [], []
    for ct in ["Basal", "Ciliated", "Secretory", "Goblet", "Fibroblast"]:
        paired = vector(x235, m235, "CRSwNP", "Nasal_polyp", "Ethmoid", ct, paired=True)
        unpaired = vector(x276, m276, "CRSwNP", "Nasal_polyp", "Inferior_turbinate", ct)
        if paired is not None:
            paired.insert(0, "celltype", ct); paired.insert(0, "contrast", "GSE235711_polyp_vs_ethmoid_paired")
            rows.append(paired)
        if unpaired is not None:
            unpaired.insert(0, "celltype", ct); unpaired.insert(0, "contrast", "GSE276503_polyp_vs_IT")
            rows.append(unpaired)
        if paired is not None and unpaired is not None:
            z = paired[["gene", "effect"]].merge(unpaired[["gene", "effect"]], on="gene", suffixes=("_paired", "_unpaired"))
            rho, p = spearmanr(z.effect_paired, z.effect_unpaired)
            strong = z[(z.effect_paired.abs() >= .5) & (z.effect_unpaired.abs() >= .5)]
            concord.append({"comparison": "polyp_transition_cross_study", "celltype": ct,
                            "genes": len(z), "rho": rho, "p": p, "strong": len(strong),
                            "same_direction_rate": (np.sign(strong.effect_paired) == np.sign(strong.effect_unpaired)).mean()})

        # Is the polyp transition merely an amplification of the mucosal disease vector?
        d = disease[(disease.dataset == "GSE235711") & (disease.case == "CRSwNP") &
                    (disease.celltype == ct)][["gene", "log2FC"]]
        if paired is not None and len(d):
            z = paired[["gene", "effect"]].merge(d, on="gene")
            rho, p = spearmanr(z.effect, z.log2FC)
            concord.append({"comparison": "paired_polyp_transition_vs_ethmoid_disease", "celltype": ct,
                            "genes": len(z), "rho": rho, "p": p, "strong": np.nan,
                            "same_direction_rate": np.nan})

    all_effects = pd.concat(rows, ignore_index=True)
    all_effects.to_csv(OUT / "anatomical_transition_gene_effects.csv.gz", index=False, compression="gzip")
    summary = pd.DataFrame(concord)
    summary.to_csv(OUT / "anatomical_transition_concordance.csv", index=False)

    core = ["EREG", "CCL20", "TNFAIP3", "FOSL1", "CXCL8", "HRH1", "CST1", "SLC26A4",
            "SCGB1A1", "FOXJ1", "PIFO", "MUC5AC", "POSTN", "COL1A1", "COL3A1", "TGFBI"]
    all_effects[all_effects.gene.isin(core)].to_csv(OUT / "core_gene_anatomical_effects.csv", index=False)

    sensitivity = []
    for threshold in [20, 50, 100]:
        for ct in ["Basal", "Ciliated", "Secretory", "Goblet", "Fibroblast"]:
            pv = vector(x235, m235, "CRSwNP", "Nasal_polyp", "Ethmoid", ct,
                        paired=True, cell_min=threshold)
            uv = vector(x276, m276, "CRSwNP", "Nasal_polyp", "Inferior_turbinate", ct,
                        cell_min=threshold)
            if pv is None or uv is None:
                continue
            z = pv[["gene", "effect"]].merge(uv[["gene", "effect"]], on="gene", suffixes=("_paired", "_unpaired"))
            rho, p = spearmanr(z.effect_paired, z.effect_unpaired)
            strong = z[(z.effect_paired.abs() >= .5) & (z.effect_unpaired.abs() >= .5)]
            sensitivity.append({"cell_min": threshold, "celltype": ct,
                                "n_paired": int(pv.n_a.iloc[0]), "n_polyp": int(uv.n_a.iloc[0]),
                                "n_it": int(uv.n_b.iloc[0]), "rho": rho, "p": p,
                                "strong": len(strong),
                                "same_direction_rate": (np.sign(strong.effect_paired) == np.sign(strong.effect_unpaired)).mean()})
    sensitivity = pd.DataFrame(sensitivity)
    sensitivity.to_csv(OUT / "cell_threshold_sensitivity.csv", index=False)
    print(summary.to_string(index=False))
    print("\nTHRESHOLD SENSITIVITY")
    print(sensitivity.to_string(index=False))
    print("\nPAIRED TOP CONSISTENT")
    q = all_effects[(all_effects.contrast.str.contains("paired")) & (all_effects.sign_consistency == 1)]
    print(q.assign(abs_effect=q.effect.abs()).sort_values(["celltype", "abs_effect"], ascending=[True, False]).groupby("celltype").head(12).to_string(index=False))


if __name__ == "__main__":
    main()
