from pathlib import Path
from io import StringIO
import csv, gzip, sys

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from validate_serous_bulk_crs import read_geo, annotation, collapse

RAW = ROOT / "data" / "raw" / "bulk_crs"
OUT = ROOT / "results" / "advanced_coexpression"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260831
RNG = np.random.default_rng(SEED)


def load_gse136825():
    counts = pd.read_csv(RAW / "GSE136825_genecounts_20190903.txt.gz", sep="\t", index_col=0)
    counts = counts[counts.columns[5:]]
    manifest = pd.read_csv(ROOT / "results" / "preflight" / "sample_manifest.csv")
    feature_path = manifest.loc[manifest.dataset == "GSE235711", "features"].iloc[0]
    feat = pd.read_csv(feature_path, sep="\t", header=None)
    mapping = pd.Series(feat.iloc[:, 1].values, index=feat.iloc[:, 0].str.split(".").str[0]).drop_duplicates()
    counts.index = counts.index.str.split(".").str[0]
    counts["symbol"] = mapping.reindex(counts.index).fillna("").values
    counts = counts[counts.symbol.ne("")].groupby("symbol").sum()
    expr = np.log2(counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1) * 1e6 + 0.5)
    _, meta = read_geo(RAW / "GSE136825_series_matrix.txt.gz")
    sources = np.array(meta["!Sample_source_name_ch1"][0])
    groups = np.array(["CRSwNP_NP" if s == "Nasal Polyp Tissue" else
                       "CRSwNP_IT" if s == "Nasal Polyp Inferior Turbinate" else "Healthy_IT" for s in sources])
    return expr, groups


def load_gse36830():
    expr, meta = read_geo(RAW / "GSE36830_series_matrix.txt.gz")
    expr = collapse(expr, annotation(RAW / "GPL570.annot.gz"))
    sources = np.array(meta["!Sample_source_name_ch1"][0])
    titles = np.array(meta["!Sample_title"][0])
    groups = []
    for source, title in zip(sources, titles):
        text = (source + " " + title).lower()
        if "control" in text and "uncinate" in text:
            groups.append("Healthy_UT")
        elif "without nasal polyps" in text or "crssnp" in text:
            groups.append("CRSsNP_UT")
        elif "nasal polyp tissue" in text or ("polyp" in text and "uncinate" not in text):
            groups.append("CRSwNP_NP")
        elif "with nasal polyps" in text or "crswnp" in text:
            groups.append("CRSwNP_UT")
        else:
            groups.append(title)
    return expr, np.array(groups)


def residualize(expr, groups):
    x = expr.copy()
    for group in np.unique(groups):
        cols = np.flatnonzero(groups == group)
        x.iloc[:, cols] = x.iloc[:, cols].sub(x.iloc[:, cols].mean(axis=1), axis=0)
    return x


def zscore_rows(x):
    return x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=1).replace(0, np.nan), axis=0)


def corr_with_score(expr, score):
    x = expr.to_numpy(float)
    y = np.asarray(score, float)
    x = x - np.nanmean(x, axis=1, keepdims=True)
    y = y - np.nanmean(y)
    den = np.sqrt(np.nansum(x * x, axis=1) * np.nansum(y * y))
    return np.divide(np.nansum(x * y, axis=1), den, out=np.full(x.shape[0], np.nan), where=den > 0)


def permutation_density(unit_rows, n_genes, n_perm=500):
    vals = []
    universe = np.arange(unit_rows.shape[0])
    for _ in range(n_perm):
        idx = RNG.choice(universe, size=n_genes, replace=False)
        u = unit_rows[idx]
        c = u @ u.T
        vals.append(np.nanmean(np.abs(c[np.triu_indices(n_genes, 1)])))
    return np.asarray(vals)


def matched_permutation_density(unit_rows, genes, module, mean_expression, residual_variance,
                                n_perm=5000, n_bins=5, seed=SEED):
    """Expression- and residual-variance-matched validation density null.

    Genes are assigned to a 5 x 5 rank-quantile grid using validation-cohort
    mean expression and within-group residual variance. Each random set exactly
    reproduces the frozen module's stratum counts, excludes module genes, and is
    sampled without replacement within each stratum.
    """
    genes = pd.Index(genes)
    metrics = pd.DataFrame({
        "mean_expression": pd.Series(mean_expression, index=genes, dtype=float),
        "residual_variance": pd.Series(residual_variance, index=genes, dtype=float),
    })
    valid = (np.isfinite(metrics.mean_expression) & np.isfinite(metrics.residual_variance) &
             (metrics.residual_variance > 0))
    metrics = metrics.loc[valid].copy()
    metrics["expression_bin"] = pd.qcut(
        metrics.mean_expression.rank(method="first"), q=n_bins, labels=False
    ).astype(int)
    metrics["variance_bin"] = pd.qcut(
        metrics.residual_variance.rank(method="first"), q=n_bins, labels=False
    ).astype(int)
    metrics["stratum"] = metrics.expression_bin.astype(str) + "_" + metrics.variance_bin.astype(str)

    missing = pd.Index(module).difference(metrics.index)
    if len(missing):
        raise ValueError(f"Module genes missing finite matched-null metrics: {missing.tolist()}")

    module_counts = metrics.loc[module, "stratum"].value_counts().sort_index()
    candidates = metrics.loc[~metrics.index.isin(module)].copy()
    pools = {}
    audit_rows = []
    for stratum, needed in module_counts.items():
        pool_genes = candidates.index[candidates.stratum == stratum]
        if len(pool_genes) < needed:
            raise ValueError(f"Matched stratum {stratum} has {len(pool_genes)} candidates for {needed} genes")
        pools[stratum] = genes.get_indexer(pool_genes)
        audit_rows.append({
            "stratum": stratum,
            "module_genes": int(needed),
            "candidate_genes": int(len(pool_genes)),
        })

    rng = np.random.default_rng(seed)
    vals = np.empty(n_perm, dtype=float)
    tri = np.triu_indices(len(module), 1)
    for b in range(n_perm):
        chosen = []
        for stratum, needed in module_counts.items():
            chosen.extend(rng.choice(pools[stratum], size=int(needed), replace=False).tolist())
        u = unit_rows[np.asarray(chosen, dtype=int)]
        c = u @ u.T
        vals[b] = np.nanmean(np.abs(c[tri]))
        if (b + 1) % 500 == 0:
            print(f"matched permutation {b + 1}/{n_perm}", flush=True)

    audit = pd.DataFrame(audit_rows)
    audit["matching_expression_bins"] = n_bins
    audit["matching_variance_bins"] = n_bins
    audit["eligible_nonmodule_universe"] = len(candidates)
    audit["seed"] = seed
    return vals, audit, metrics


def main():
    d_expr, d_group = load_gse136825()
    v_expr, v_group = load_gse36830()
    common = d_expr.index.intersection(v_expr.index)
    d_expr, v_expr = d_expr.loc[common], v_expr.loc[common]
    locked = pd.read_csv(ROOT / "results" / "locked_40_gene_module.csv").gene
    locked = [g for g in locked if g in common]

    d_res, v_res = residualize(d_expr, d_group), residualize(v_expr, v_group)
    d_locked_score = zscore_rows(d_res.loc[locked]).mean(axis=0)
    v_locked_score = zscore_rows(v_res.loc[locked]).mean(axis=0)
    d_cor = pd.Series(corr_with_score(d_res, d_locked_score), index=common, name="discovery_module_cor")
    v_cor = pd.Series(corr_with_score(v_res, v_locked_score), index=common, name="validation_module_cor")

    # Freeze the positive discovery neighborhood, excluding the 40 locked genes.
    neighborhood = d_cor.drop(index=locked).dropna().sort_values(ascending=False).head(160).index.tolist()
    module = locked + neighborhood
    node = pd.concat([d_cor.loc[module], v_cor.loc[module]], axis=1).reset_index(names="gene")
    node["is_locked"] = node.gene.isin(locked)
    node["same_positive_direction"] = (node.discovery_module_cor > 0) & (node.validation_module_cor > 0)

    dc = np.corrcoef(d_res.loc[module].to_numpy())
    vc = np.corrcoef(v_res.loc[module].to_numpy())
    tri = np.triu_indices(len(module), 1)
    edge_preservation_rho, edge_preservation_p = spearmanr(dc[tri], vc[tri], nan_policy="omit")
    discovery_density = float(np.nanmean(np.abs(dc[tri])))
    validation_density = float(np.nanmean(np.abs(vc[tri])))

    # Validation density is compared with random same-sized gene sets.
    v_arr = v_res.to_numpy(float)
    v_arr = v_arr - np.nanmean(v_arr, axis=1, keepdims=True)
    v_norm = np.sqrt(np.nansum(v_arr * v_arr, axis=1, keepdims=True))
    v_unit = np.divide(v_arr, v_norm, out=np.zeros_like(v_arr), where=v_norm > 0)
    perm = permutation_density(v_unit, len(module), 500)
    empirical_p = (1 + np.sum(perm >= validation_density)) / (1 + len(perm))

    # Stronger sensitivity null: preserve the module's joint distribution of
    # validation mean expression and within-group residual variance.
    matched_perm, matched_audit, matched_metrics = matched_permutation_density(
        v_unit,
        common,
        module,
        mean_expression=v_expr.mean(axis=1).reindex(common).to_numpy(float),
        residual_variance=v_res.var(axis=1, ddof=1).reindex(common).to_numpy(float),
        n_perm=5000,
        n_bins=5,
        seed=SEED,
    )
    matched_empirical_p = (1 + np.sum(matched_perm >= validation_density)) / (1 + len(matched_perm))

    # Retain the complete 198-choose-2 edge universe separately from the
    # thresholded consensus network. Cross-cohort edge correlation is computed
    # over all 19,503 pairs; only is_consensus pairs are used to build the
    # readable positive co-expression network.
    all_edges = pd.DataFrame({
        "source": np.asarray(module, dtype=object)[tri[0]],
        "target": np.asarray(module, dtype=object)[tri[1]],
        "cor_discovery": dc[tri],
        "cor_validation": vc[tri],
    })
    all_edges["is_consensus"] = (
        (all_edges.cor_discovery >= 0.55) &
        (all_edges.cor_validation >= 0.30)
    )
    all_edges["consensus_weight"] = np.nan
    consensus_mask = all_edges.is_consensus
    all_edges.loc[consensus_mask, "consensus_weight"] = np.sqrt(
        all_edges.loc[consensus_mask, "cor_discovery"] *
        all_edges.loc[consensus_mask, "cor_validation"]
    )
    edges = all_edges.loc[
        all_edges.is_consensus,
        ["source", "target", "cor_discovery", "cor_validation", "consensus_weight"],
    ].copy()
    if len(edges):
        degree = pd.concat([edges.source, edges.target]).value_counts()
        weight = pd.concat([edges[["source", "consensus_weight"]].rename(columns={"source": "gene"}),
                            edges[["target", "consensus_weight"]].rename(columns={"target": "gene"})])
        strength = weight.groupby("gene").consensus_weight.sum()
        node["consensus_degree"] = node.gene.map(degree).fillna(0).astype(int)
        node["consensus_strength"] = node.gene.map(strength).fillna(0.0)
    else:
        node["consensus_degree"] = 0
        node["consensus_strength"] = 0.0

    summary = pd.DataFrame([{
        "discovery_dataset": "GSE136825", "validation_dataset": "GSE36830",
        "locked_genes": len(locked), "discovery_neighbors": len(neighborhood),
        "module_genes": len(module), "validation_same_positive_fraction": node.same_positive_direction.mean(),
        "gene_module_correlation_rho": spearmanr(node.discovery_module_cor, node.validation_module_cor).statistic,
        "edge_correlation_preservation_rho": edge_preservation_rho,
        "edge_correlation_preservation_p": edge_preservation_p,
        "discovery_absolute_density": discovery_density, "validation_absolute_density": validation_density,
        "validation_density_empirical_p": empirical_p,
        "validation_density_matched_empirical_p": matched_empirical_p,
        "matched_permutations": len(matched_perm),
        "matched_expression_bins": 5,
        "matched_variance_bins": 5,
        "matched_nonmodule_universe": int(matched_audit.eligible_nonmodule_universe.iloc[0]),
        "permutation_seed": SEED,
        "consensus_edges": len(edges)
    }])
    node.sort_values(["consensus_degree", "consensus_strength"], ascending=False).to_csv(OUT / "preserved_module_nodes.csv", index=False)
    edges.sort_values("consensus_weight", ascending=False).to_csv(OUT / "preserved_module_edges.csv", index=False)
    all_edges.sort_values(["is_consensus", "cor_discovery", "cor_validation"],
                          ascending=[False, False, False]).to_csv(
        OUT / "all_module_gene_pairs.csv", index=False
    )
    summary.to_csv(OUT / "module_preservation_summary.csv", index=False)
    pd.DataFrame({"permuted_validation_density": perm}).to_csv(OUT / "module_density_permutation.csv", index=False)
    pd.DataFrame({"matched_permuted_validation_density": matched_perm}).to_csv(
        OUT / "module_density_matched_permutation.csv", index=False
    )
    matched_audit.to_csv(OUT / "module_density_matched_strata_audit.csv", index=False)
    matched_metrics.reset_index(names="gene").to_csv(
        OUT / "module_density_matching_metrics.csv", index=False
    )
    print(summary.to_string(index=False))
    print("\nTop preserved hubs")
    print(node.sort_values(["consensus_degree", "consensus_strength"], ascending=False).head(30).to_string(index=False))


if __name__ == "__main__":
    main()
