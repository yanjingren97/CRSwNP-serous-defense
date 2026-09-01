from pathlib import Path
from collections import defaultdict
import sys

ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse
from celltype_preflight import genes_from, BROAD

MANIFEST = ROOT / "results" / "preflight" / "sample_manifest.csv"
OUT = ROOT / "results" / "secretory_subtypes"
OUT.mkdir(parents=True, exist_ok=True)

# Reference-guided panels fixed before looking at disease effects.
PANELS = {
    "Surface_club": ["SCGB1A1", "SCGB3A1", "KRT4", "KRT13", "KLF3", "CYP2F1"],
    "Serous_glandular": ["LTF", "LYZ", "BPIFB1", "SLPI", "WFDC2", "PRR4", "LPO"],
    "Mucous_glandular": ["MUC5B", "BPIFB2", "AGR2", "SPDEF", "ZG16B", "PIP"],
    "Goblet": ["MUC5AC", "CLCA1", "SPDEF", "AGR2", "FCGBP"],
    "Inflammatory_secretory": ["CST1", "SLC26A4", "SERPINB3", "SERPINB4", "CCL26", "POSTN"],
    "Transitional_secretory": ["KRT4", "KRT13", "KRT17", "KRT19", "SFN", "KRT16"],
}
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT4", "KRT13"]
SECRETORY = sorted(set(sum(PANELS.values(), [])))
IMMUNE = ["PTPRC", "LST1", "TYROBP", "CD3D", "CD79A", "NKG7"]
IG_PREFIX = ("IGH", "IGK", "IGL")
MARKERS = sorted(set(EPI + SECRETORY + IMMUNE + sum(BROAD.values(), [])))


def rows_for(genes, wanted):
    return {g: np.flatnonzero(genes == g) for g in wanted if np.any(genes == g)}


def one(row):
    genes = genes_from(row.features)
    x = mmread(row.matrix).tocsr()
    if x.shape[0] != len(genes) and x.shape[1] == len(genes):
        x = x.T.tocsr()
    total = np.asarray(x.sum(axis=0)).ravel()
    ngene = np.asarray((x > 0).sum(axis=0)).ravel()
    mtidx = np.flatnonzero(np.char.startswith(genes.astype(str), "MT-"))
    mt = np.asarray(x[mtidx].sum(axis=0)).ravel() if len(mtidx) else np.zeros_like(total)
    keep = (ngene >= 200) & (total > 0) & (100 * mt / np.maximum(total, 1) < 20)
    x, total = x[:, keep], total[keep]

    idx = rows_for(genes, MARKERS)
    vals = {}
    det = {}
    for g, ii in idx.items():
        raw = np.asarray(x[ii].sum(axis=0)).ravel()
        vals[g] = np.log1p(raw * 1e4 / total)
        det[g] = raw > 0

    def score(panel):
        present = [g for g in panel if g in vals]
        if not present:
            return np.zeros(x.shape[1]), np.zeros(x.shape[1], int)
        return np.mean([vals[g] for g in present], axis=0), np.sum([det[g] for g in present], axis=0)

    sec_s, sec_d = score(SECRETORY)
    broad_scores, broad_detected = [], []
    for panel in BROAD.values():
        a, b = score(panel); broad_scores.append(a); broad_detected.append(b)
    broad_scores = np.vstack(broad_scores).T
    broad_detected = np.vstack(broad_detected).T
    broad_best = np.argmax(broad_scores, axis=1)
    epithelial_index = list(BROAD).index("Epithelial")
    epithelial_wins = broad_best == epithelial_index
    epithelial_score = broad_scores[:, epithelial_index]
    epithelial_detected = broad_detected[:, epithelial_index]
    # A cell must be epithelial by the same four-compartment competition used
    # in the broad audit; a lone ambient secretory transcript is insufficient.
    candidates = epithelial_wins & (epithelial_detected >= 1) & (epithelial_score >= .15) & (sec_d >= 1)

    names = list(PANELS)
    ss, dd = [], []
    for panel in PANELS.values():
        a, b = score(panel); ss.append(a); dd.append(b)
    ss, dd = np.vstack(ss).T, np.vstack(dd).T
    best = np.argmax(ss, axis=1)
    ordered = np.sort(ss, axis=1)
    margin = ordered[:, -1] - ordered[:, -2]
    labels = np.full(x.shape[1], "Not_secretory", dtype=object)
    accepted = candidates & (dd[np.arange(len(best)), best] >= 1) & (ss[np.arange(len(best)), best] >= .12)
    labels[accepted] = np.asarray(names, object)[best[accepted]]
    ambiguous = accepted & (margin < .05)
    labels[ambiguous] = "Ambiguous_secretory"

    # Ambient-RNA flags are sample-level diagnostics, not automatic exclusions.
    igidx = np.flatnonzero(np.array([str(g).startswith(IG_PREFIX) for g in genes]))
    ig = np.asarray(x[igidx].sum(axis=0)).ravel() if len(igidx) else np.zeros(x.shape[1])
    composition, sums = [], {}
    codes, unique = pd.factorize(genes, sort=False)
    mapper = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(unique), len(codes)))
    for lab in names + ["Ambiguous_secretory"]:
        mask = labels == lab
        composition.append({"dataset": row.dataset, "sample": row.sample, "donor": row.donor,
                            "disease": row.disease, "tissue": row.tissue, "subtype": lab,
                            "n_cells": int(mask.sum()), "n_secretory_candidates": int(candidates.sum()),
                            "candidate_fraction": float(mask.sum()/max(candidates.sum(), 1)),
                            "median_ig_fraction": float(np.median(ig[mask]/np.maximum(total[mask], 1))) if mask.any() else np.nan})
        rawsum = np.asarray(x[:, mask].sum(axis=1)).ravel().astype(np.int64)
        sums[lab] = (np.asarray(mapper @ rawsum).ravel().astype(np.int64), int(mask.sum()))
    return np.asarray(unique), composition, sums, int(candidates.sum())


def main():
    manifest = pd.read_csv(MANIFEST)
    composition, errors = [], []
    for dataset, rows in manifest.groupby("dataset", sort=False):
        print("DATASET", dataset, flush=True)
        universe = np.asarray(sorted(set.intersection(*[set(genes_from(p)) for p in rows.features])))
        donor_counts, donor_cells, donor_meta = {}, defaultdict(int), {}
        for row in rows.itertuples(index=False):
            print(" ", row.sample, flush=True)
            try:
                genes, comp, sums, nc = one(row)
                composition.extend(comp)
                order = pd.Index(genes).get_indexer(universe)
                for subtype, (counts, cells) in sums.items():
                    key = (row.donor, row.disease, row.tissue, subtype)
                    donor_counts[key] = donor_counts.get(key, 0) + counts[order]
                    donor_cells[key] += cells
                    donor_meta[key] = {"dataset": dataset, "donor": row.donor, "disease": row.disease,
                                       "tissue": row.tissue, "subtype": subtype}
            except Exception as exc:
                errors.append({"dataset": dataset, "sample": row.sample, "error": f"{type(exc).__name__}: {exc}"})
        keys = list(donor_counts)
        mat = np.column_stack([donor_counts[k] for k in keys])
        cols = ["|".join(k) for k in keys]
        pd.DataFrame(mat, index=universe, columns=cols).to_csv(OUT/f"{dataset}_subtype_counts.tsv.gz", sep="\t", compression="gzip")
        pd.DataFrame([{**donor_meta[k], "column": "|".join(k), "n_cells": donor_cells[k]} for k in keys]).to_csv(OUT/f"{dataset}_subtype_metadata.csv", index=False)
    pd.DataFrame(composition).to_csv(OUT/"sample_subtype_composition.csv", index=False)
    pd.DataFrame(errors).to_csv(OUT/"errors.csv", index=False)
    d = pd.DataFrame(composition).groupby(["dataset", "donor", "disease", "tissue", "subtype"], as_index=False).agg(
        n_cells=("n_cells", "sum"), n_secretory_candidates=("n_secretory_candidates", "sum"),
        median_ig_fraction=("median_ig_fraction", "median"))
    d["n_classified_secretory"] = d.groupby(["dataset", "donor", "disease", "tissue"]).n_cells.transform("sum")
    d["proportion_of_candidates"] = d.n_cells / d.n_secretory_candidates.clip(lower=1)
    d["proportion_of_classified"] = d.n_cells / d.n_classified_secretory.clip(lower=1)
    # Backward-compatible alias; all new formal models use the explicitly named columns.
    d["proportion"] = d["proportion_of_classified"]
    d.to_csv(OUT/"donor_subtype_composition.csv", index=False)
    print(d.groupby(["dataset", "disease", "tissue", "subtype"]).agg(donors=("donor", "nunique"), median_cells=("n_cells", "median"), median_prop=("proportion", "median")).reset_index().to_string(index=False))


if __name__ == "__main__":
    main()
