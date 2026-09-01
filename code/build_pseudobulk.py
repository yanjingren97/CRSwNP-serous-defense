from pathlib import Path
from collections import defaultdict
import sys

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse

from celltype_preflight import BROAD, EPI, ALL_MARKERS, genes_from, marker_matrix, scores

MANIFEST = ROOT / "results" / "preflight" / "sample_manifest.csv"
OUT = ROOT / "results" / "pseudobulk"
OUT.mkdir(parents=True, exist_ok=True)
TARGET_TYPES = ["Basal", "Ciliated", "Secretory", "Goblet", "Fibroblast"]


def classify_and_sum(row):
    genes = genes_from(row.features)
    x = mmread(row.matrix).tocsr()
    if x.shape[0] != len(genes) and x.shape[1] == len(genes):
        x = x.T.tocsr()
    total = np.asarray(x.sum(axis=0)).ravel()
    n_gene = np.asarray((x > 0).sum(axis=0)).ravel()
    mt_idx = np.flatnonzero(np.char.startswith(genes.astype(str), "MT-"))
    mt = np.asarray(x[mt_idx].sum(axis=0)).ravel() if len(mt_idx) else np.zeros_like(total)
    keep = (n_gene >= 200) & (total > 0) & (100 * mt / np.maximum(total, 1) < 20)
    x = x[:, keep]
    total = total[keep]
    rawm, marker_names = marker_matrix(x, genes)
    logm = np.log1p(rawm * (1e4 / total)[None, :])

    bs, bd = scores(logm, marker_names, BROAD)
    bi = np.argmax(bs, axis=1)
    broad = np.asarray(list(BROAD), dtype=object)[bi]
    broad[(bs[np.arange(len(bi)), bi] < .15) | (bd[np.arange(len(bi)), bi] == 0)] = "Unresolved"

    fine = np.full(x.shape[1], "Not_epithelial", dtype=object)
    emask = broad == "Epithelial"
    if emask.any():
        es, ed = scores(logm[:, emask], marker_names, EPI)
        ei = np.argmax(es, axis=1)
        lab = np.asarray(list(EPI), dtype=object)[ei]
        lab[(es[np.arange(len(ei)), ei] < .12) | (ed[np.arange(len(ei)), ei] == 0)] = "Epithelial_other"
        fine[emask] = lab

    ans = {}
    for ct in TARGET_TYPES:
        mask = (broad == "Fibroblast") if ct == "Fibroblast" else (fine == ct)
        ans[ct] = (np.asarray(x[:, mask].sum(axis=1)).ravel().astype(np.int64), int(mask.sum()))
    # Collapse duplicate symbols within each sample before cross-sample alignment.
    codes, unique = pd.factorize(genes, sort=False)
    mapper = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                               shape=(len(unique), len(codes)))
    ans = {ct: (np.asarray(mapper @ val).ravel().astype(np.int64), nc)
           for ct, (val, nc) in ans.items()}
    return np.asarray(unique), ans


def collapse_duplicate_genes(genes, matrix):
    # Preserve first-seen order while summing duplicated gene symbols.
    codes, unique = pd.factorize(genes, sort=False)
    mapper = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                               shape=(len(unique), len(codes)))
    return np.asarray(unique), mapper @ matrix


def main():
    manifest = pd.read_csv(MANIFEST)
    if len(sys.argv) > 1:
        manifest = manifest[manifest.dataset.isin(sys.argv[1:])]
    for dataset, rows in manifest.groupby("dataset", sort=False):
        print(f"DATASET {dataset}", flush=True)
        gene_sets = [set(genes_from(p)) for p in rows.features]
        base_genes = np.asarray(sorted(set.intersection(*gene_sets)))
        print(f"  common gene universe: {len(base_genes)}", flush=True)
        donor_counts = {}
        donor_cells = defaultdict(int)
        metadata = {}
        errors = []
        for row in rows.itertuples(index=False):
            print(f"  {row.sample}", flush=True)
            try:
                genes, vals = classify_and_sum(row)
                order = pd.Index(genes).get_indexer(base_genes)
                if np.any(order < 0):
                    raise ValueError("Common-gene alignment failed")
                vals = {ct: (count[order], nc) for ct, (count, nc) in vals.items()}
                for ct, (count, nc) in vals.items():
                    key = (row.donor, row.disease, row.tissue, ct)
                    donor_counts[key] = donor_counts.get(key, 0) + count
                    donor_cells[key] += nc
                    metadata[key] = {"dataset": dataset, "donor": row.donor, "disease": row.disease,
                                     "tissue": row.tissue, "celltype": ct}
            except Exception as exc:
                errors.append({"dataset": dataset, "sample": row.sample,
                               "error": f"{type(exc).__name__}: {exc}"})
                print(f"    SKIP {errors[-1]['error']}", flush=True)
        keys = list(donor_counts)
        mat = np.column_stack([donor_counts[k] for k in keys])
        cols = ["|".join(k) for k in keys]
        pd.DataFrame(mat, index=base_genes, columns=cols).to_csv(
            OUT / f"{dataset}_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")
        meta = pd.DataFrame([{**metadata[k], "column": "|".join(k), "n_cells": donor_cells[k]} for k in keys])
        meta.to_csv(OUT / f"{dataset}_pseudobulk_metadata.csv", index=False)
        if errors:
            pd.DataFrame(errors).to_csv(OUT / f"{dataset}_errors.csv", index=False)
        print(f"  wrote {mat.shape[0]} genes x {mat.shape[1]} pseudobulks", flush=True)


if __name__ == "__main__":
    main()
