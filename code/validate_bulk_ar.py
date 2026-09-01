from pathlib import Path
import csv
import gzip
import io
import sys

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

RAW = ROOT / "data" / "raw" / "bulk_ar"
OUT = ROOT / "results" / "preflight"

SIGNATURES = {
    "type2_epithelial": ["POSTN", "CCL26", "CST1", "CLCA1", "ALOX15", "IL13RA2", "SERPINB2"],
    "epithelial_remodeling": ["KRT17", "KRT6A", "KRT16", "MMP9", "TGFBI", "LAMC2", "ITGA6"],
    "ecm_fibrotic": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "SPARC"],
    "interferon": ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "STAT1"],
}

CONFIG = {
    "GSE19187": {"gpl": "GPL6244", "select": lambda x: x.str.startswith(("Rhinitis", "Healthy")),
                 "group": lambda x: np.where(x.str.startswith("Rhinitis"), "AR", "Healthy")},
    "GSE43523": {"gpl": "GPL6883", "select": lambda x: pd.Series(True, index=x.index),
                 "group": lambda x: np.where(x.str.contains("rhinitis", case=False), "AR", "Healthy")},
    "GSE44037": {"gpl": "GPL13158", "select": lambda x: x.str.contains("nasal_epithelium") & x.str.contains("rhinitis|healthy"),
                 "group": lambda x: np.where(x.str.contains("rhinitis"), "AR", "Healthy")},
}


def parse_annot(gpl):
    path = RAW / f"{gpl}.annot.gz"
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                break
        tab = pd.read_csv(f, sep="\t", dtype=str)
    tab = tab[["ID", "Gene symbol"]].dropna()
    tab["Gene symbol"] = tab["Gene symbol"].str.split("///").str[0].str.strip().str.upper()
    return tab[tab["Gene symbol"] != ""]


def parse_series(gse):
    path = RAW / f"{gse}_series_matrix.txt.gz"
    titles = None
    table_lines = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                table_lines.append(line)
    expr = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    title = pd.Series(titles, index=expr.columns, name="title")
    return expr.apply(pd.to_numeric, errors="coerce"), title


def hedges_g(case, control):
    n1, n0 = len(case), len(control)
    pooled = np.sqrt(((n1 - 1) * np.var(case, ddof=1) + (n0 - 1) * np.var(control, ddof=1)) / (n1 + n0 - 2))
    if not np.isfinite(pooled) or pooled == 0:
        return np.nan
    d = (np.mean(case) - np.mean(control)) / pooled
    return d * (1 - 3 / (4 * (n1 + n0) - 9))


def main():
    all_scores, results = [], []
    for gse, cfg in CONFIG.items():
        expr, title = parse_series(gse)
        keep = cfg["select"](title)
        expr, title = expr.loc[:, keep], title[keep]
        groups = pd.Series(cfg["group"](title), index=title.index)

        annot = parse_annot(cfg["gpl"])
        merged = annot.merge(expr, left_on="ID", right_index=True, how="inner")
        gene_expr = merged.drop(columns="ID").groupby("Gene symbol").mean(numeric_only=True)

        score_df = pd.DataFrame(index=expr.columns)
        score_df["dataset"] = gse
        score_df["sample"] = score_df.index
        score_df["title"] = title
        score_df["group"] = groups
        for sig, genes in SIGNATURES.items():
            present = [g for g in genes if g in gene_expr.index]
            vals = gene_expr.loc[present].T
            z = (vals - vals.mean()) / vals.std(ddof=1).replace(0, np.nan)
            score_df[sig] = z.mean(axis=1)
            case = score_df.loc[score_df.group == "AR", sig].dropna().to_numpy()
            ctrl = score_df.loc[score_df.group == "Healthy", sig].dropna().to_numpy()
            p = mannwhitneyu(case, ctrl, alternative="two-sided", method="auto").pvalue
            results.append({"dataset": gse, "signature": sig, "n_AR": len(case), "n_Healthy": len(ctrl),
                            "genes_present": len(present), "genes_total": len(genes),
                            "mean_difference_z": np.mean(case) - np.mean(ctrl),
                            "median_difference_z": np.median(case) - np.median(ctrl),
                            "hedges_g": hedges_g(case, ctrl), "mannwhitney_p": p})
        all_scores.append(score_df.reset_index(drop=True))

    scores = pd.concat(all_scores, ignore_index=True)
    res = pd.DataFrame(results)
    scores.to_csv(OUT / "bulk_ar_signature_scores.csv", index=False)
    res.to_csv(OUT / "bulk_ar_signature_validation.csv", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
