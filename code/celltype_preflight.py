from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
from scipy.io import mmread

INP = ROOT / "results" / "preflight" / "sample_manifest.csv"
OUT = ROOT / "results" / "celltype_preflight"
OUT.mkdir(parents=True, exist_ok=True)

BROAD = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT4", "KRT5"],
    "Immune": ["PTPRC", "LST1", "TYROBP", "CD3D", "CD79A", "MS4A1", "NKG7"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "CLDN5", "RAMP2"],
}

EPI = {
    "Basal": ["KRT5", "KRT14", "KRT15", "TP63", "KRT17"],
    "Ciliated": ["FOXJ1", "PIFO", "TPPP3", "CAPS", "CETN2"],
    "Secretory": ["SCGB1A1", "SCGB3A1", "KRT4", "KRT13", "BPIFB1"],
    "Goblet": ["MUC5AC", "SPDEF", "AGR2", "CLCA1"],
    "Ionocyte": ["FOXI1", "CFTR", "ASCL3"],
    "Tuft": ["POU2F3", "TRPM5", "AVIL"],
    "Cycling": ["MKI67", "TOP2A", "UBE2C"],
}

SIGNATURES = {
    "type2_epithelial": ["POSTN", "CCL26", "CST1", "CLCA1", "ALOX15", "IL13RA2", "SERPINB2"],
    "epithelial_remodeling": ["KRT17", "KRT6A", "KRT16", "MMP9", "TGFBI", "LAMC2", "ITGA6"],
    "ecm_fibrotic": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "SPARC"],
}
ALL_MARKERS = sorted(set(sum(BROAD.values(), []) + sum(EPI.values(), []) + sum(SIGNATURES.values(), [])))


def genes_from(path):
    x = pd.read_csv(path, sep="\t", header=None)
    return (x.iloc[:, 1] if x.shape[1] > 1 else x.iloc[:, 0]).astype(str).str.upper().to_numpy()


def marker_matrix(x, genes):
    rows, names = [], []
    for g in ALL_MARKERS:
        idx = np.flatnonzero(genes == g)
        if len(idx):
            rows.append(x[idx].sum(axis=0))
            names.append(g)
    if not rows:
        raise ValueError("No marker genes found")
    return np.vstack([np.asarray(v).ravel() for v in rows]), np.asarray(names)


def scores(logm, names, panels):
    ans = np.zeros((logm.shape[1], len(panels)), dtype=np.float32)
    detected = np.zeros_like(ans, dtype=np.int16)
    for j, genes in enumerate(panels.values()):
        idx = np.flatnonzero(np.isin(names, genes))
        if len(idx):
            ans[:, j] = logm[idx].mean(axis=0)
            detected[:, j] = (logm[idx] > 0).sum(axis=0)
    return ans, detected


def one(row):
    genes = genes_from(row.features)
    x = mmread(row.matrix).tocsr()
    if x.shape[0] != len(genes) and x.shape[1] == len(genes):
        x = x.T.tocsr()
    total = np.asarray(x.sum(axis=0)).ravel()
    n_gene = np.asarray((x > 0).sum(axis=0)).ravel()
    mt_idx = np.flatnonzero(np.char.startswith(genes.astype(str), "MT-"))
    mt = np.asarray(x[mt_idx].sum(axis=0)).ravel() if len(mt_idx) else np.zeros_like(total)
    keep = (n_gene >= 200) & (total > 0) & (100 * mt / np.maximum(total, 1) < 20)
    if keep.sum() < 50:
        raise ValueError(f"Only {keep.sum()} cells pass QC")
    x = x[:, keep]
    total = total[keep]
    rawm, names = marker_matrix(x, genes)
    logm = np.log1p(rawm * (1e4 / total)[None, :])

    bscore, bdet = scores(logm, names, BROAD)
    bi = np.argmax(bscore, axis=1)
    broad_names = np.asarray(list(BROAD))
    broad = broad_names[bi]
    # Require at least one detected marker and a modest normalized score.
    mx = bscore[np.arange(len(bi)), bi]
    nd = bdet[np.arange(len(bi)), bi]
    broad[(mx < 0.15) | (nd == 0)] = "Unresolved"

    epi_mask = broad == "Epithelial"
    epi_labels = np.full(len(broad), "Not_epithelial", dtype=object)
    if epi_mask.any():
        escore, edet = scores(logm[:, epi_mask], names, EPI)
        ei = np.argmax(escore, axis=1)
        labels = np.asarray(list(EPI), dtype=object)[ei]
        emx = escore[np.arange(len(ei)), ei]
        end = edet[np.arange(len(ei)), ei]
        labels[(emx < 0.12) | (end == 0)] = "Epithelial_other"
        epi_labels[epi_mask] = labels

    comp = []
    for level, labels in [("broad", broad), ("epithelial", epi_labels[epi_mask])]:
        vals, cnt = np.unique(labels, return_counts=True)
        denom = len(labels)
        for lab, n in zip(vals, cnt):
            comp.append({"dataset": row.dataset, "sample": row.sample, "donor": row.donor,
                         "disease": row.disease, "tissue": row.tissue, "level": level,
                         "celltype": lab, "n_cells": int(n), "proportion": n / denom})

    pseudo = []
    for compartment, mask in [("Epithelial", epi_mask), ("Fibroblast", broad == "Fibroblast")]:
        for sig, glist in SIGNATURES.items():
            idx = np.flatnonzero(np.isin(names, glist))
            counts = rawm[idx][:, mask].sum(axis=1) if len(idx) and mask.any() else np.zeros(len(idx))
            pseudo.append({"dataset": row.dataset, "sample": row.sample, "donor": row.donor,
                           "disease": row.disease, "tissue": row.tissue, "compartment": compartment,
                           "signature": sig, "n_cells": int(mask.sum()),
                           "marker_counts": float(np.sum(counts)),
                           "all_counts": float(total[mask].sum()) if mask.any() else 0.0,
                           "genes_present": int(len(idx))})
    return comp, pseudo


def main():
    manifest = pd.read_csv(INP)
    comp, pseudo, errors = [], [], []
    for i, row in enumerate(manifest.itertuples(index=False), 1):
        print(f"[{i}/{len(manifest)}] {row.dataset} {row.sample}", flush=True)
        try:
            a, b = one(row)
            comp.extend(a); pseudo.extend(b)
        except Exception as exc:
            errors.append({"dataset": row.dataset, "sample": row.sample,
                           "error": f"{type(exc).__name__}: {exc}"})
            print(f"  SKIP {errors[-1]['error']}", flush=True)
    comp = pd.DataFrame(comp)
    pseudo = pd.DataFrame(pseudo)
    comp.to_csv(OUT / "sample_celltype_composition.csv", index=False)
    pseudo.to_csv(OUT / "sample_compartment_signatures.csv", index=False)
    pd.DataFrame(errors).to_csv(OUT / "errors.csv", index=False)

    # Donor/tissue summary: technical libraries aggregate by cell counts.
    d = comp.groupby(["dataset", "donor", "disease", "tissue", "level", "celltype"], as_index=False).n_cells.sum()
    d["proportion"] = d.n_cells / d.groupby(["dataset", "donor", "disease", "tissue", "level"]).n_cells.transform("sum")
    d.to_csv(OUT / "donor_celltype_composition.csv", index=False)
    summary = d.groupby(["dataset", "disease", "tissue", "level", "celltype"]).agg(
        donors=("donor", "nunique"), median_cells=("n_cells", "median"),
        median_proportion=("proportion", "median")).reset_index()
    summary.to_csv(OUT / "cohort_celltype_summary.csv", index=False)
    print(summary[summary.level == "broad"].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
