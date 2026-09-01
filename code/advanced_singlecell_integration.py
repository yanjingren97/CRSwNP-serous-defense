from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.sparse.linalg import eigsh
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.manifold import TSNE
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, f1_score
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
import anndata as ad
import harmonypy as hm

SEED = 20260831
RNG = np.random.default_rng(SEED)
OUT = ROOT / "results" / "advanced_singlecell"
OUT.mkdir(parents=True, exist_ok=True)

BROAD = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT4", "KRT5"],
    "Immune": ["PTPRC", "LST1", "TYROBP", "CD3D", "CD79A", "MS4A1", "NKG7"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "CLDN5", "RAMP2"],
}
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
DEFINITION = set(sum(BROAD.values(), []) + sum(PANELS.values(), []) + EPI)
CILIA = ["FOXJ1", "PIFO", "TPPP3", "CAPS", "CETN2", "DNAH5", "DNAI2", "RSPH1", "CFAP54", "IFT88"]


def genes_from(path):
    x = pd.read_csv(path, sep="\t", header=None)
    return (x.iloc[:, 1] if x.shape[1] > 1 else x.iloc[:, 0]).astype(str).str.upper().to_numpy()


def choose_feature_panel():
    mats = []
    for ds in ["GSE235711", "GSE276503"]:
        x = pd.read_csv(ROOT / f"results/secretory_subtypes/{ds}_subtype_counts.tsv.gz", sep="\t", index_col=0)
        lib = x.sum(axis=0).replace(0, np.nan)
        mats.append(np.log2(x.div(lib, axis=1) * 1e6 + 0.5))
    common = mats[0].index.intersection(mats[1].index)
    z = pd.concat([m.loc[common] for m in mats], axis=1)
    var = z.var(axis=1)
    eligible = [g for g in common if g not in DEFINITION and not g.startswith(("MT-", "RPL", "RPS", "IG"))]
    top = var.loc[eligible].sort_values(ascending=False).head(5000).index.tolist()
    locked = pd.read_csv(ROOT / "results/locked_40_gene_module.csv").gene.astype(str).str.upper().tolist()
    required = sorted(set(sum(BROAD.values(), []) + sum(PANELS.values(), []) + CILIA + locked))
    panel = sorted(set(top + [g for g in required if g in common]))
    pd.DataFrame({"gene": panel, "role": ["required" if g in required else "pseudobulk_variable" for g in panel]}).to_csv(OUT / "integration_gene_panel.csv", index=False)
    return np.asarray(panel)


def collapse_to_panel(x, genes, panel):
    pos = {g: i for i, g in enumerate(panel)}
    src, dst = [], []
    for i, g in enumerate(genes):
        j = pos.get(str(g))
        if j is not None:
            src.append(i); dst.append(j)
    mapper = sparse.csr_matrix((np.ones(len(src), dtype=np.float32), (dst, src)), shape=(len(panel), x.shape[0]))
    return mapper @ x


def score_panels(x, genes, panels):
    total = np.asarray(x.sum(axis=0)).ravel()
    vals, det = [], []
    index = pd.Index(genes)
    for glist in panels.values():
        ii = index.get_indexer([g for g in glist if g in index])
        ii = ii[ii >= 0]
        if len(ii):
            raw = x[ii]
            norm = raw.multiply(1e4 / np.maximum(total, 1))
            vals.append(np.asarray(np.log1p(norm).mean(axis=0)).ravel())
            det.append(np.asarray((raw > 0).sum(axis=0)).ravel())
        else:
            vals.append(np.zeros(x.shape[1])); det.append(np.zeros(x.shape[1]))
    return np.column_stack(vals), np.column_stack(det)


def read_one(row, panel):
    genes = genes_from(row.features)
    x = mmread(row.matrix).tocsr()
    if x.shape[0] != len(genes) and x.shape[1] == len(genes):
        x = x.T.tocsr()
    total = np.asarray(x.sum(axis=0)).ravel()
    ngene = np.asarray((x > 0).sum(axis=0)).ravel()
    mtidx = np.flatnonzero(np.char.startswith(genes.astype(str), "MT-"))
    mt = np.asarray(x[mtidx].sum(axis=0)).ravel() if len(mtidx) else np.zeros_like(total)
    keep = (ngene >= 200) & (total > 0) & (100 * mt / np.maximum(total, 1) < 20)
    x = x[:, keep]
    bar = pd.read_csv(row.barcodes, sep="\t", header=None).iloc[:, 0].astype(str).to_numpy()[keep]
    y = collapse_to_panel(x, genes, panel)

    bscore, bdet = score_panels(y, panel, BROAD)
    bbest = np.argmax(bscore, axis=1)
    epi_i = list(BROAD).index("Epithelial")
    epithelial = (bbest == epi_i) & (bdet[:, epi_i] >= 1) & (bscore[:, epi_i] >= 0.15)
    sscore, sdet = score_panels(y, panel, PANELS)
    pbest = np.argmax(sscore, axis=1)
    ordered = np.sort(sscore, axis=1)
    margin = ordered[:, -1] - ordered[:, -2]
    accepted = epithelial & (sdet[np.arange(len(pbest)), pbest] >= 1) & (sscore[np.arange(len(pbest)), pbest] >= 0.12)
    labels = np.full(y.shape[1], "Not_secretory", dtype=object)
    labels[accepted] = np.asarray(list(PANELS), dtype=object)[pbest[accepted]]
    labels[accepted & (margin < 0.05)] = "Ambiguous_secretory"
    candidates = accepted
    if candidates.sum() == 0:
        return None

    selected = []
    candidate_idx = np.flatnonzero(candidates)
    for lab in sorted(np.unique(labels[candidates])):
        ii = candidate_idx[labels[candidates] == lab]
        if len(ii) > 350:
            ii = RNG.choice(ii, 350, replace=False)
        selected.extend(ii.tolist())
    selected = np.asarray(sorted(selected))
    obs = pd.DataFrame({
        "cell_id": [f"{row.dataset}|{row.sample}|{b}" for b in bar[selected]],
        "dataset": row.dataset,
        "sample": row.sample,
        "donor": row.donor,
        "disease": row.disease,
        "tissue": row.tissue,
        "subtype": labels[selected],
        "total_counts": total[keep][selected],
        "n_genes": ngene[keep][selected],
        "subtype_margin": margin[selected],
    }).set_index("cell_id")
    return ad.AnnData(X=y[:, selected].T.tocsr(), obs=obs, var=pd.DataFrame(index=panel))


def mean_score(adata, genes):
    ii = adata.var_names.get_indexer([g for g in genes if g in adata.var_names])
    ii = ii[ii >= 0]
    if not len(ii):
        return np.zeros(adata.n_obs)
    return np.asarray(adata.X[:, ii].mean(axis=1)).ravel()


def sparse_log_normalize(x):
    x = x.tocsr().astype(np.float32)
    lib = np.asarray(x.sum(axis=1)).ravel()
    x = sparse.diags(1e4 / np.maximum(lib, 1)) @ x
    x.data = np.log1p(x.data)
    return x.tocsr()


def sparse_variance(x):
    mean = np.asarray(x.mean(axis=0)).ravel()
    sq = x.copy(); sq.data **= 2
    return np.asarray(sq.mean(axis=0)).ravel() - mean ** 2


def harmony(pcs, obs, key="dataset"):
    z = hm.run_harmony(pcs, obs[[key]].copy(), [key], random_state=SEED, verbose=False)
    return np.asarray(z.Z_corr).T


def integration_and_mapping(adata):
    adata.layers["counts"] = adata.X.copy()
    adata.X = sparse_log_normalize(adata.X)
    eligible = np.array([g not in DEFINITION and not g.startswith(("MT-", "RPL", "RPS", "IG")) for g in adata.var_names])
    adata.var["eligible_embedding"] = eligible
    var = sparse_variance(adata.X)
    hvg_idx = np.flatnonzero(eligible)[np.argsort(var[eligible])[-min(2500, int(eligible.sum())):]]
    adata.var["highly_variable"] = False
    adata.var.iloc[hvg_idx, adata.var.columns.get_loc("highly_variable")] = True
    pcs = TruncatedSVD(n_components=40, random_state=SEED).fit_transform(adata.X[:, hvg_idx])
    corrected = harmony(pcs, adata.obs, "dataset")
    adata.obsm["X_pca"] = pcs
    adata.obsm["X_pca_harmony"] = corrected
    # Barnes-Hut t-SNE supplies a non-linear display without altering patient-level inference.
    emb = TSNE(n_components=2, perplexity=50, learning_rate="auto", init="pca", max_iter=1200,
               method="barnes_hut", random_state=SEED, n_jobs=-1).fit_transform(corrected[:, :30])
    adata.obsm["X_tsne_harmony"] = emb
    km = KMeans(n_clusters=8, n_init=30, random_state=SEED).fit(adata.obsm["X_pca_harmony"][:, :30])
    adata.obs["unsupervised_state"] = [f"State_{x+1}" for x in km.labels_]

    rows, cms = [], []
    for train, test in [("GSE235711", "GSE276503"), ("GSE276503", "GSE235711")]:
        a = adata.obs.dataset == train; b = adata.obs.dataset == test
        clf = KNeighborsClassifier(n_neighbors=25, weights="distance")
        clf.fit(adata.obsm["X_pca_harmony"][a, :30], adata.obs.loc[a, "subtype"])
        pred = clf.predict(adata.obsm["X_pca_harmony"][b, :30])
        truth = adata.obs.loc[b, "subtype"].to_numpy()
        rows.append({"train_dataset": train, "test_dataset": test, "n_test": int(b.sum()),
                     "balanced_accuracy": balanced_accuracy_score(truth, pred),
                     "macro_f1": f1_score(truth, pred, average="macro")})
        labels = sorted(set(truth) | set(pred))
        cm = confusion_matrix(truth, pred, labels=labels, normalize="true")
        for i, t in enumerate(labels):
            for j, p in enumerate(labels):
                cms.append({"train_dataset": train, "test_dataset": test, "truth": t, "predicted": p, "fraction": cm[i, j]})
    pd.DataFrame(rows).to_csv(OUT / "cross_dataset_label_transfer_metrics.csv", index=False)
    pd.DataFrame(cms).to_csv(OUT / "cross_dataset_label_transfer_confusion.csv", index=False)


def serous_continuum(adata):
    s = adata[adata.obs.subtype == "Serous_glandular"].copy()
    # Recompute a serous-specific non-definition diffusion geometry.
    eligible = s.var.eligible_embedding.to_numpy(bool)
    var = sparse_variance(s.X)
    idx = np.flatnonzero(eligible)[np.argsort(var[eligible])[-min(1500, int(eligible.sum())):]]
    pcs = TruncatedSVD(n_components=30, random_state=SEED).fit_transform(s.X[:, idx])
    corrected = harmony(pcs, s.obs, "dataset")
    nn = NearestNeighbors(n_neighbors=25, metric="euclidean", n_jobs=-1).fit(corrected[:, :25])
    dist, ind = nn.kneighbors()
    sigma = float(np.median(dist[:, -1]))
    rows = np.repeat(np.arange(s.n_obs), ind.shape[1])
    w = np.exp(-np.square(dist.ravel()) / max(sigma ** 2, 1e-8))
    graph = sparse.csr_matrix((w, (rows, ind.ravel())), shape=(s.n_obs, s.n_obs))
    graph = graph.maximum(graph.T)
    degree = np.asarray(graph.sum(axis=1)).ravel()
    norm = sparse.diags(1 / np.sqrt(np.maximum(degree, 1e-8)))
    kernel = norm @ graph @ norm
    # Use an explicitly seeded ARPACK start vector so the diffusion eigensolver
    # is reproducible independently of NumPy's process-global random state.
    eig_v0 = np.random.default_rng(SEED).normal(size=kernel.shape[0])
    evals, evecs = eigsh(
        kernel,
        k=11,
        which="LA",
        v0=eig_v0,
        tol=1e-8,
        maxiter=max(1000, kernel.shape[0] * 5),
    )
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order][1:11], evecs[:, order][:, 1:11]
    diff = evecs * evals[None, :]
    locked = pd.read_csv(ROOT / "results/locked_40_gene_module.csv").gene.astype(str).str.upper().tolist()
    s.obs["mature_score"] = mean_score(s, locked)
    s.obs["cilia_score"] = mean_score(s, CILIA)
    ref = (s.obs.dataset == "GSE235711") & (s.obs.disease == "Healthy")
    orient = s.obs.mature_score - s.obs.cilia_score
    pool = np.flatnonzero(ref.to_numpy()) if ref.any() else np.arange(s.n_obs)
    root = int(pool[np.argmax(orient.iloc[pool].to_numpy())])
    d = np.sqrt(np.square(diff - diff[root]).sum(axis=1))
    s.obs["failure_axis"] = pd.Series(d).rank(method="average", pct=True).to_numpy()
    cols = ["dataset", "sample", "donor", "disease", "tissue", "subtype", "mature_score", "cilia_score", "failure_axis"]
    s.obs[cols].to_csv(OUT / "serous_cell_state_axis.csv")
    donor = s.obs.reset_index().groupby(["dataset", "donor", "disease", "tissue"], as_index=False).agg(
        n_cells=("cell_id", "size"), median_failure_axis=("failure_axis", "median"),
        median_mature_score=("mature_score", "median"), median_cilia_score=("cilia_score", "median"))
    donor.to_csv(OUT / "serous_donor_state_axis.csv", index=False)
    # Save diffusion coordinates for plotting and cross-dataset geometry checks.
    diff_df = pd.DataFrame(diff[:, :5], index=s.obs_names, columns=[f"DC{i}" for i in range(1, 6)])
    pd.concat([s.obs[cols], diff_df], axis=1).to_csv(OUT / "serous_diffusion_coordinates.csv")


def main():
    panel = choose_feature_panel()
    manifest = pd.read_csv(ROOT / "results/preflight/sample_manifest.csv")
    manifest = manifest[manifest.dataset.isin(["GSE235711", "GSE276503"])]
    objects, errors, audit = [], [], []
    for i, row in enumerate(manifest.itertuples(index=False), 1):
        print(f"[{i}/{len(manifest)}] {row.dataset} {row.sample}", flush=True)
        try:
            obj = read_one(row, panel)
            if obj is not None:
                objects.append(obj)
                audit.append({"dataset": row.dataset, "sample": row.sample, "cells_retained_balanced": obj.n_obs})
        except Exception as exc:
            errors.append({"dataset": row.dataset, "sample": row.sample, "error": f"{type(exc).__name__}: {exc}"})
            print(" SKIP", errors[-1]["error"], flush=True)
    pd.DataFrame(errors).to_csv(OUT / "integration_errors.csv", index=False)
    pd.DataFrame(audit).to_csv(OUT / "integration_sample_audit.csv", index=False)
    if not objects:
        raise RuntimeError("No cell objects were constructed")
    adata = ad.concat(objects, join="inner", merge="same", index_unique=None)
    adata.obs_names_make_unique()
    integration_and_mapping(adata)
    coords = pd.DataFrame(adata.obsm["X_tsne_harmony"], index=adata.obs_names, columns=["tSNE1", "tSNE2"])
    pd.concat([adata.obs, coords], axis=1).to_csv(OUT / "integrated_secretory_embedding.csv")
    serous_continuum(adata)
    # Pandas 3 may use Arrow-backed strings, which h5py cannot serialize directly.
    adata.obs.index = pd.Index(adata.obs.index.astype(str).tolist(), dtype=object, name="cell_id")
    for col in adata.obs.columns:
        if pd.api.types.is_string_dtype(adata.obs[col].dtype):
            adata.obs[col] = adata.obs[col].astype(str).astype("category")
    adata.write_h5ad(OUT / "integrated_secretory_balanced.h5ad", compression="gzip")
    print("FINAL", adata.shape, flush=True)


if __name__ == "__main__":
    main()
