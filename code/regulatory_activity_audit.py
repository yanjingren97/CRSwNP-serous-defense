from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp, ttest_ind, norm, chi2

from secretory_subtype_pseudobulk_validation import load

RES = ROOT / "data" / "resources" / "regulatory"
OUT = ROOT / "results" / "advanced_regulatory"
OUT.mkdir(parents=True, exist_ok=True)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(np.nan_to_num(p, nan=1.0))
    q = np.full(n, np.nan)
    last = 1.0
    for rank, idx in reversed(list(enumerate(order, 1))):
        if np.isnan(p[idx]):
            continue
        last = min(last, p[idx] * n / rank)
        q[idx] = last
    return q


def read_progeny():
    raw = pd.read_csv(RES / "PROGENy_annotations.tsv", sep="\t", dtype=str)
    x = raw.pivot_table(index=["genesymbol", "record_id"], columns="label", values="value", aggfunc="first").reset_index()
    x["weight"] = pd.to_numeric(x["weight"], errors="coerce")
    x["p_value"] = pd.to_numeric(x["p_value"], errors="coerce")
    x = x.dropna(subset=["pathway", "weight", "p_value"])
    x = x[x.p_value < 0.05].rename(columns={"pathway": "source", "genesymbol": "target"})
    # Keep the strongest 500 responsive genes per pathway to limit diffuse signatures.
    x = x.sort_values("p_value").groupby("source", observed=True).head(500)
    return x[["source", "target", "weight", "p_value"]].drop_duplicates(["source", "target"])


def read_collectri():
    x = pd.read_csv(RES / "CollecTRI_regulons.csv")
    x = x[["source", "target", "weight", "resources", "references"]].dropna(subset=["source", "target", "weight"])
    x["weight"] = pd.to_numeric(x.weight, errors="coerce")
    return x.dropna(subset=["weight"]).drop_duplicates(["source", "target"])


def activity_matrix(expr, network, min_targets=10):
    """Weighted mean of gene-wise standardized logCPM, normalized by L2 weight.

    This is a transparent signed-regulon enrichment score. Gene standardization is
    performed within dataset and the same frozen network is used in both cohorts.
    """
    common = sorted(set(expr.index).intersection(network.target))
    z = expr.loc[common].copy()
    mu = z.mean(axis=1)
    sd = z.std(axis=1, ddof=1).replace(0, np.nan)
    z = z.sub(mu, axis=0).div(sd, axis=0).fillna(0.0)
    acts, audit = {}, []
    for src, g in network[network.target.isin(common)].groupby("source", observed=True):
        g = g.drop_duplicates("target")
        if len(g) < min_targets:
            continue
        w = g.set_index("target").weight.astype(float)
        w = w.loc[w.index.intersection(z.index)]
        score = z.loc[w.index].T.dot(w) / np.sqrt(np.square(w).sum())
        acts[src] = score
        audit.append({"source": src, "n_targets_measured": len(w), "weight_l2": np.sqrt(np.square(w).sum())})
    return pd.DataFrame(acts).T, pd.DataFrame(audit)


def select_serous(expr, meta):
    m = meta[(meta.subtype == "Serous_glandular") & (meta.n_cells >= 20)].copy()
    return expr[m.index], m


def contrast_paired(a, m):
    z = m[(m.disease == "CRSwNP") & m.tissue.isin(["Nasal_polyp", "Ethmoid"])].copy()
    z["column_id"] = z.index
    w = z.reset_index(drop=True).pivot(index="donor", columns="tissue", values="column_id").dropna()
    np_cols, et_cols = w.Nasal_polyp.tolist(), w.Ethmoid.tolist()
    d = a[np_cols].to_numpy() - a[et_cols].to_numpy()
    eff = d.mean(axis=1)
    sd = d.std(axis=1, ddof=1)
    se = sd / np.sqrt(d.shape[1])
    p = ttest_1samp(d, 0, axis=1, nan_policy="omit").pvalue
    # Fraction of leave-one-pair-out mean effects retaining the full-data sign.
    loo = []
    for i in range(d.shape[1]):
        loo.append(np.delete(d, i, axis=1).mean(axis=1))
    loo = np.asarray(loo).T
    stability = (np.sign(loo) == np.sign(eff[:, None])).mean(axis=1)
    return pd.DataFrame({"source": a.index, "effect_paired": eff, "se_paired": se,
                         "p_paired": p, "n_paired": d.shape[1], "loo_sign_stability": stability})


def contrast_unpaired(a, m):
    z = m[(m.disease == "CRSwNP") & m.tissue.isin(["Nasal_polyp", "Inferior_turbinate"])]
    np_cols = z[z.tissue == "Nasal_polyp"].index.tolist()
    it_cols = z[z.tissue == "Inferior_turbinate"].index.tolist()
    y1, y0 = a[np_cols].to_numpy(), a[it_cols].to_numpy()
    eff = y1.mean(axis=1) - y0.mean(axis=1)
    v1, v0 = y1.var(axis=1, ddof=1) / len(np_cols), y0.var(axis=1, ddof=1) / len(it_cols)
    se = np.sqrt(v1 + v0)
    p = ttest_ind(y1, y0, axis=1, equal_var=False, nan_policy="omit").pvalue
    return pd.DataFrame({"source": a.index, "effect_external": eff, "se_external": se,
                         "p_external": p, "n_external_np": len(np_cols), "n_external_it": len(it_cols)})


def synthesize(paired, external, kind):
    z = paired.merge(external, on="source")
    w1 = 1 / z.se_paired.pow(2).replace(0, np.nan)
    w2 = 1 / z.se_external.pow(2).replace(0, np.nan)
    z["meta_effect"] = (w1 * z.effect_paired + w2 * z.effect_external) / (w1 + w2)
    z["meta_se"] = np.sqrt(1 / (w1 + w2))
    z["meta_p"] = 2 * norm.sf(np.abs(z.meta_effect / z.meta_se))
    z["meta_fdr"] = bh(z.meta_p)
    z["direction_concordant"] = np.sign(z.effect_paired) == np.sign(z.effect_external)
    z["fisher_p"] = chi2.sf(-2 * (np.log(z.p_paired.clip(lower=1e-300)) + np.log(z.p_external.clip(lower=1e-300))), 4)
    z["fisher_fdr"] = bh(z.fisher_p)
    z.insert(0, "network", kind)
    return z.sort_values(["direction_concordant", "fisher_fdr"], ascending=[False, True])


def run_one(name, net, min_targets):
    x1, m1 = load("GSE235711")
    x2, m2 = load("GSE276503")
    x1, m1 = select_serous(x1, m1)
    x2, m2 = select_serous(x2, m2)
    a1, q1 = activity_matrix(x1, net, min_targets=min_targets)
    a2, q2 = activity_matrix(x2, net, min_targets=min_targets)
    common = a1.index.intersection(a2.index)
    a1, a2 = a1.loc[common], a2.loc[common]
    q = q1.merge(q2, on="source", suffixes=("_GSE235711", "_GSE276503"))
    q.insert(0, "network", name)
    q.to_csv(OUT / f"{name.lower()}_target_coverage.csv", index=False)

    s1 = a1.T.reset_index(names="column").melt(id_vars="column", var_name="source", value_name="activity")
    s1 = s1.merge(m1.reset_index(names="column"), on="column"); s1.insert(0, "network", name)
    s2 = a2.T.reset_index(names="column").melt(id_vars="column", var_name="source", value_name="activity")
    s2 = s2.merge(m2.reset_index(names="column"), on="column"); s2.insert(0, "network", name)
    pd.concat([s1, s2], ignore_index=True).to_csv(OUT / f"{name.lower()}_sample_activities.csv", index=False)

    result = synthesize(contrast_paired(a1, m1), contrast_unpaired(a2, m2), name)
    if name == "CollecTRI":
        detection_cut = np.log2(1.0 + 0.5)  # CPM >= 1 in the log2(CPM + 0.5) matrix
        s1 = pd.Series({g: float((x1.loc[g] >= detection_cut).mean()) if g in x1.index else 0.0
                        for g in result.source}, name="source_detected_fraction_paired")
        s2 = pd.Series({g: float((x2.loc[g] >= detection_cut).mean()) if g in x2.index else 0.0
                        for g in result.source}, name="source_detected_fraction_external")
        result = result.merge(s1.rename_axis("source").reset_index(), on="source")
        result = result.merge(s2.rename_axis("source").reset_index(), on="source")
    else:
        result["source_detected_fraction_paired"] = np.nan
        result["source_detected_fraction_external"] = np.nan
    result.to_csv(OUT / f"{name.lower()}_cross_cohort_contrasts.csv", index=False)
    return result


def main():
    ct = read_collectri()
    pg = read_progeny()
    ct.to_csv(OUT / "collectri_network_frozen.csv.gz", index=False, compression="gzip")
    pg.to_csv(OUT / "progeny_network_frozen.csv.gz", index=False, compression="gzip")
    tf = run_one("CollecTRI", ct, 20)
    pw = run_one("PROGENy", pg, 25)
    combined = pd.concat([tf, pw], ignore_index=True)
    combined["replicated_priority"] = (combined.direction_concordant &
                                        (combined.loo_sign_stability >= 0.8) &
                                        (combined.p_paired < 0.10) &
                                        (combined.p_external < 0.05) &
                                        (combined.fisher_fdr < 0.10))
    is_tf = combined.network == "CollecTRI"
    expression_supported = ((combined.source_detected_fraction_paired >= 0.25) &
                            (combined.source_detected_fraction_external >= 0.25))
    combined.loc[is_tf & ~expression_supported, "replicated_priority"] = False
    combined.to_csv(OUT / "regulatory_cross_cohort_summary.csv", index=False)
    print("\nReplicated regulatory candidates")
    cols = ["network", "source", "effect_paired", "effect_external", "p_paired",
            "p_external", "fisher_fdr", "loo_sign_stability", "meta_effect"]
    print(combined[combined.replicated_priority][cols].sort_values("fisher_fdr").to_string(index=False))
    print("\nAudit counts")
    print(combined.groupby("network").agg(n_tested=("source", "size"),
          n_concordant=("direction_concordant", "sum"), n_priority=("replicated_priority", "sum")).to_string())


if __name__ == "__main__":
    main()
