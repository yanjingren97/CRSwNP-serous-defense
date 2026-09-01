from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
from scipy.io import mmread

MANIFEST = ROOT / "results" / "preflight" / "sample_manifest.csv"
OUT = ROOT / "results" / "preflight"

SIGNATURES = {
    "type2_epithelial": ["POSTN", "CCL26", "CST1", "CLCA1", "ALOX15", "IL13RA2", "SERPINB2"],
    "epithelial_remodeling": ["KRT17", "KRT6A", "KRT16", "MMP9", "TGFBI", "LAMC2", "ITGA6"],
    "ecm_fibrotic": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "SPARC"],
    "interferon": ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "STAT1"],
}
TARGETS = sorted(set(sum(SIGNATURES.values(), [])))


def read_genes(path):
    tab = pd.read_csv(path, sep="\t", header=None)
    symbols = tab.iloc[:, 1] if tab.shape[1] > 1 else tab.iloc[:, 0]
    return symbols.astype(str).str.upper().to_numpy()


def sample_counts(row):
    genes = read_genes(row.features)
    x = mmread(row.matrix).tocsr()
    if x.shape[0] != len(genes) and x.shape[1] == len(genes):
        x = x.T.tocsr()
    if x.shape[0] != len(genes):
        raise ValueError(f"shape {x.shape} incompatible with {len(genes)} genes")
    totals = np.asarray(x.sum(axis=1)).ravel()
    lib = float(totals.sum())
    result = {g: float(totals[genes == g].sum()) for g in TARGETS}
    result.update({"library_size": lib, "n_cells": x.shape[1]})
    return result


def main():
    manifest = pd.read_csv(MANIFEST)
    rows = []
    for i, row in enumerate(manifest.itertuples(index=False), 1):
        print(f"[{i}/{len(manifest)}] {row.sample}", flush=True)
        try:
            vals = sample_counts(row)
            vals.update({k: getattr(row, k) for k in ["dataset", "sample", "donor", "disease", "tissue"]})
            rows.append(vals)
        except Exception as exc:
            print(f"  SKIP {type(exc).__name__}: {exc}", flush=True)
    sample = pd.DataFrame(rows)
    sample.to_csv(OUT / "target_gene_sample_counts.csv", index=False)

    # Aggregate technical libraries and paired tissues only at the explicit
    # donor/disease/tissue level; do not invent pairing for GSE276503.
    keys = ["dataset", "donor", "disease", "tissue"]
    donor = sample.groupby(keys, as_index=False)[TARGETS + ["library_size", "n_cells"]].sum()
    for gene in TARGETS:
        donor[gene] = np.log1p(1e6 * donor[gene] / donor.library_size)
    for name, genes in SIGNATURES.items():
        donor[name] = donor[genes].mean(axis=1)
    donor.to_csv(OUT / "donor_signature_scores.csv", index=False)

    contrasts = [
        ("GSE235711", "Ethmoid", "CRSsNP", "Healthy"),
        ("GSE235711", "Ethmoid", "CRSwNP", "Healthy"),
        ("GSE235711", "Nasal_polyp", "CRSwNP", None),
        ("GSE276503", "Inferior_turbinate", "CRSwNP", "Healthy"),
        ("GSE276503", "Nasal_polyp", "CRSwNP", None),
        ("GSE261706", "Inferior_turbinate", "AR", "Healthy"),
    ]
    out = []
    for ds, tissue, case, control in contrasts:
        a = donor[(donor.dataset == ds) & (donor.tissue == tissue) & (donor.disease == case)]
        b = donor[(donor.dataset == ds) & (donor.tissue == tissue) & (donor.disease == control)] if control else None
        for sig in SIGNATURES:
            item = {"dataset": ds, "tissue": tissue, "case": case, "control": control or "descriptive_only",
                    "signature": sig, "n_case": len(a), "case_median": a[sig].median()}
            if b is not None:
                item.update({"n_control": len(b), "control_median": b[sig].median(),
                             "median_difference": a[sig].median() - b[sig].median()})
            out.append(item)
    result = pd.DataFrame(out)
    result.to_csv(OUT / "signature_contrasts.csv", index=False)
    print(result.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
