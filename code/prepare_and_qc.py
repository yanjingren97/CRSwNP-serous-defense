import csv
import gzip
import re
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
from scipy.io import mmread

EXTRACTED = ROOT / "data" / "derived" / "extracted"
PREPARED = ROOT / "data" / "derived" / "prepared"
RESULTS = ROOT / "results" / "preflight"
PREPARED.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)


def unpack_gse235711():
    source = EXTRACTED / "GSE235711"
    target = PREPARED / "GSE235711"
    target.mkdir(exist_ok=True)
    for archive in sorted(source.glob("*.tar.gz")):
        sample = archive.name.replace("_filtered.tar.gz", "")
        sample_dir = target / sample
        matrix = sample_dir / "matrix.mtx.gz"
        if matrix.exists():
            continue
        sample_dir.mkdir(exist_ok=True)
        with tarfile.open(archive, "r:gz") as tar:
            safe_members = []
            for member in tar.getmembers():
                if Path(member.name).name in {"matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"}:
                    member.name = Path(member.name).name
                    safe_members.append(member)
            tar.extractall(sample_dir, members=safe_members, filter="data")


def manifest_rows():
    rows = []
    # GSE235711: sample titles encode group, site and donor number.
    for d in sorted((PREPARED / "GSE235711").iterdir()):
        name = re.sub(r"^GSM\d+_", "", d.name)
        if name.startswith("Control_"):
            disease, tissue = "Healthy", "Ethmoid"
            donor = "HC_" + name.split("_")[1]
        elif name.startswith("CRSsNP_Eth_"):
            disease, tissue = "CRSsNP", "Ethmoid"
            donor = "CRSsNP_" + name.split("_")[-1]
        elif name.startswith("CRSwNP_Eth_"):
            disease, tissue = "CRSwNP", "Ethmoid"
            donor = "CRSwNP_" + name.split("_")[-1]
        elif name.startswith("CRSwNP_NP_"):
            disease, tissue = "CRSwNP", "Nasal_polyp"
            donor = "CRSwNP_" + name.split("_")[-1]
        else:
            raise ValueError(name)
        rows.append(dict(dataset="GSE235711", sample=d.name, donor=donor,
                         disease=disease, tissue=tissue, library="single",
                         matrix=str(d / "matrix.mtx.gz"),
                         features=str(d / "features.tsv.gz"),
                         barcodes=str(d / "barcodes.tsv.gz")))

    # GSE276503: publication reports 15 CRSwNP and 2 HC donors, but GEO does not
    # expose a cross-tissue donor key. Keep every biopsy separate and never infer pairing.
    d = EXTRACTED / "GSE276503"
    for matrix in sorted(d.glob("*_matrix.mtx.gz")):
        prefix = matrix.name.replace("_matrix.mtx.gz", "")
        prefix = prefix.replace("GSM8499604", "GSM8499604_CRSNP33") if prefix == "GSM8499604" else prefix
        label = re.sub(r"^GSM\d+_", "", prefix)
        if label.startswith("HCIT"):
            disease, tissue = "Healthy", "Inferior_turbinate"
        elif label.startswith("CRSIT"):
            disease, tissue = "CRSwNP", "Inferior_turbinate"
        elif label.startswith("CRSMT"):
            disease, tissue = "CRSwNP", "Middle_turbinate"
        elif label.startswith("CRSNP"):
            disease, tissue = "CRSwNP", "Nasal_polyp"
        else:
            raise ValueError(label)
        base = matrix.name.replace("_matrix.mtx.gz", "")
        genes = d / f"{base}_genes.tsv.gz"
        barcodes = d / f"{base}_barcodes.tsv.gz"
        if not genes.exists() and "GSM8499604" in base:
            genes = d / "GSM8499604_CRSNP33_genes.tsv.gz"
            barcodes = d / "GSM8499604_CRSNP33_barcodes.tsv.gz"
        rows.append(dict(dataset="GSE276503", sample=label, donor=label,
                         disease=disease, tissue=tissue, library="single",
                         matrix=str(matrix), features=str(genes), barcodes=str(barcodes)))

    # GSE261706: two lanes per donor; donor, not lane, is the biological replicate.
    d = EXTRACTED / "GSE261706"
    for matrix in sorted(d.glob("*_matrix.mtx.gz")):
        base = matrix.name.replace("_filtered_matrix.mtx.gz", "")
        gsm, library_id = base.split("_", 1)
        # mapping from GEO sample order to biological donor
        gsm_num = int(gsm.replace("GSM", ""))
        donor_map = {
            8149310: ("AR1", "AR"), 8149311: ("AR1", "AR"),
            8149312: ("AR2", "AR"), 8149313: ("AR2", "AR"),
            8149314: ("CON1", "Healthy"), 8149315: ("CON1", "Healthy"),
            8149316: ("CON2", "Healthy"), 8149317: ("CON2", "Healthy"),
        }
        donor, disease = donor_map[gsm_num]
        rows.append(dict(dataset="GSE261706", sample=gsm, donor=donor,
                         disease=disease, tissue="Inferior_turbinate", library=library_id,
                         matrix=str(matrix),
                         features=str(d / f"{base}_filtered_features.tsv.gz"),
                         barcodes=str(d / f"{base}_filtered_barcodes.tsv.gz")))
    return pd.DataFrame(rows)


def read_genes(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    return np.array([r[1] if len(r) > 1 else r[0] for r in rows], dtype=object)


def count_lines_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        return sum(1 for _ in fh)


def qc_one(row):
    x = mmread(row.matrix).tocsr()  # genes x cells
    genes = read_genes(row.features)
    n_barcodes = count_lines_gz(row.barcodes)
    if x.shape[0] != len(genes) or x.shape[1] != n_barcodes:
        raise ValueError(f"Dimension mismatch {row.sample}: {x.shape}, {len(genes)}, {n_barcodes}")
    totals = np.asarray(x.sum(axis=0)).ravel()
    detected = np.diff(x.tocsc().indptr)
    mt = np.char.startswith(genes.astype(str), "MT-")
    mt_counts = np.asarray(x[mt].sum(axis=0)).ravel() if mt.any() else np.zeros(x.shape[1])
    pct_mt = np.divide(mt_counts * 100.0, totals, out=np.zeros_like(mt_counts, dtype=float), where=totals > 0)
    return dict(dataset=row.dataset, sample=row.sample, donor=row.donor,
                disease=row.disease, tissue=row.tissue, library=row.library,
                n_genes=x.shape[0], n_cells=x.shape[1], total_umis=int(totals.sum()),
                median_umis=float(np.median(totals)), median_genes=float(np.median(detected)),
                median_pct_mt=float(np.median(pct_mt)),
                cells_ge_200_genes=int((detected >= 200).sum()),
                cells_lt_20pct_mt=int((pct_mt < 20).sum()))


def main():
    unpack_gse235711()
    manifest = manifest_rows()
    manifest.to_csv(RESULTS / "sample_manifest.csv", index=False)
    qc = []
    for i, row in enumerate(manifest.itertuples(index=False), 1):
        print(f"[{i}/{len(manifest)}] {row.dataset} {row.sample}", flush=True)
        try:
            qc.append(qc_one(row))
        except Exception as exc:
            # A malformed author-supplied matrix should be documented and
            # excluded, rather than aborting the audit of every other sample.
            qc.append({
                "dataset": row.dataset,
                "sample": row.sample,
                "donor": row.donor,
                "disease": row.disease,
                "tissue": row.tissue,
                "qc_error": f"{type(exc).__name__}: {exc}",
            })
            print(f"  QC_ERROR: {type(exc).__name__}: {exc}", flush=True)
    qc = pd.DataFrame(qc)
    qc.to_csv(RESULTS / "sample_qc.csv", index=False)
    summary = (qc.groupby(["dataset", "disease", "tissue"], dropna=False)
                 .agg(samples=("sample", "count"), donors=("donor", "nunique"),
                      cells=("n_cells", "sum"), cells_ge_200_genes=("cells_ge_200_genes", "sum"),
                      median_sample_genes=("median_genes", "median"),
                      median_sample_pct_mt=("median_pct_mt", "median"))
                 .reset_index())
    summary.to_csv(RESULTS / "cohort_qc_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
