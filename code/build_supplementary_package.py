from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D


OUT = ROOT / "supplementary"
FIG = OUT / "figures"
BUILD = ROOT / "supplementary" / "build"
SEED = 20260831

NAVY = "#17324D"
BLUE = "#2F6F9F"
TEAL = "#2A9D8F"
ORANGE = "#E67E22"
RED = "#C94C4C"
GREY = "#6B7280"
LIGHT = "#EEF3F7"


BROAD = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT4", "KRT5"],
    "Immune": ["PTPRC", "LST1", "TYROBP", "CD3D", "CD79A", "MS4A1", "NKG7"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "CLDN5", "RAMP2"],
}
SECRETORY = {
    "Surface_club": ["SCGB1A1", "SCGB3A1", "KRT4", "KRT13", "KLF3", "CYP2F1"],
    "Serous_glandular": ["LTF", "LYZ", "BPIFB1", "SLPI", "WFDC2", "PRR4", "LPO"],
    "Mucous_glandular": ["MUC5B", "BPIFB2", "AGR2", "SPDEF", "ZG16B", "PIP"],
    "Goblet": ["MUC5AC", "CLCA1", "SPDEF", "AGR2", "FCGBP"],
    "Inflammatory_secretory": ["CST1", "SLC26A4", "SERPINB3", "SERPINB4", "CCL26", "POSTN"],
    "Transitional_secretory": ["KRT4", "KRT13", "KRT17", "KRT19", "SFN", "KRT16"],
}
CILIA = ["FOXJ1", "PIFO", "TPPP3", "CAPS", "CETN2", "DNAH5", "DNAI2", "RSPH1", "CFAP54", "IFT88"]


def read(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path)


def clean_value(x):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def frame_payload(df: pd.DataFrame, description: str) -> dict:
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    rows = [[clean_value(v) for v in row] for row in df.itertuples(index=False, name=None)]
    return {"description": description, "columns": list(df.columns), "rows": rows}


def union_frames(parts: list[tuple[str, pd.DataFrame]]) -> pd.DataFrame:
    out = []
    for kind, df in parts:
        z = df.copy()
        z.insert(0, "record_type", kind)
        out.append(z)
    return pd.concat(out, ignore_index=True, sort=False)


def gse136825_pairing_audit() -> pd.DataFrame:
    """Reconstruct the exact title-derived candidate-pair audit used in R."""
    scores = read("results/bulk_crs_validation/GSE136825_scores.csv")
    scores = scores[scores.group.isin(["CRSwNP_NP", "CRSwNP_IT"])].copy()
    prefix_map = {
        "PY": ("polyp", "PY-SCP"),
        "SCP": ("inferior_turbinate", "PY-SCP"),
        "TT": ("polyp", "TT-SCT"),
        "SCT": ("inferior_turbinate", "TT-SCT"),
    }
    scores["title_prefix"] = scores.title.astype(str).str.extract(r"^([^ ]+)")[0]
    scores["pair_role"] = scores.title_prefix.map(lambda x: prefix_map.get(x, (None, None))[0])
    scores["pair_series"] = scores.title_prefix.map(lambda x: prefix_map.get(x, (None, None))[1])
    scores["pair_id"] = scores.apply(
        lambda r: str(r.title)[len(str(r.title_prefix)) + 1:] if pd.notna(r.pair_role) else None,
        axis=1,
    )
    scores["pair_key"] = scores.pair_series.astype(str) + "|" + scores.pair_id.astype(str)
    scores["duplicate_within_role_key"] = scores.duplicated(
        ["pair_role", "pair_series", "pair_id"], keep=False
    )
    counterpart = scores[["pair_role", "pair_series", "pair_id", "sample", "title"]].copy()
    counterpart["pair_role"] = counterpart.pair_role.map({
        "polyp": "inferior_turbinate", "inferior_turbinate": "polyp"
    })
    counterpart = counterpart.rename(columns={
        "sample": "counterpart_sample", "title": "counterpart_title"
    })
    scores = scores.merge(
        counterpart,
        on=["pair_role", "pair_series", "pair_id"],
        how="left",
        validate="one_to_one",
    )
    scores["pair_included"] = scores.counterpart_sample.notna()
    scores["pair_status"] = np.where(scores.pair_included, "matched", "unmatched_excluded")
    scores["pairing_source_field"] = "GEO !Sample_title"
    scores["duplicate_handling"] = (
        "No duplicate role/series/identifier keys were observed; no duplicate-resolution rule was invoked"
    )
    columns = [
        "sample", "title", "source", "group", "score", "title_prefix", "pair_role",
        "pair_series", "pair_id", "pair_key", "counterpart_sample", "counterpart_title",
        "pair_included", "pair_status", "duplicate_within_role_key", "pairing_source_field",
        "duplicate_handling",
    ]
    return scores[columns].sort_values(["pair_series", "pair_id", "pair_role"], na_position="last")


def marker_table() -> pd.DataFrame:
    rows = []
    for panel, genes in BROAD.items():
        rows.extend({"entry_type": "marker", "marker_family": "broad_compartment", "panel": panel,
                     "gene": g, "predefined_before_disease_testing": True} for g in genes)
    for panel, genes in SECRETORY.items():
        rows.extend({"entry_type": "marker", "marker_family": "secretory_subtype", "panel": panel,
                     "gene": g, "predefined_before_disease_testing": True} for g in genes)
    rows.extend({"entry_type": "marker", "marker_family": "cilia_state", "panel": "Cilia", "gene": g,
                 "predefined_before_disease_testing": True} for g in CILIA)
    rules = [
        ("cell_qc", "detected_genes", ">=", 200, "Applied to every cell before annotation"),
        ("cell_qc", "mitochondrial_fraction_percent", "<", 20, "Cells with zero library size were removed"),
        ("broad_annotation", "winning_panel_mean_log_normalized_score", ">=", 0.15, "At least one panel marker detected"),
        ("secretory_annotation", "winning_subtype_mean_log_normalized_score", ">=", 0.12, "At least one subtype marker detected"),
        ("secretory_annotation", "best_minus_second_best_score_margin", ">=", 0.05, "Otherwise Ambiguous_secretory"),
        ("integration_sampling", "maximum_cells_per_sample_by_subtype", "<=", 350, "Without replacement; seed 20260831"),
        ("pseudobulk_primary", "cells_per_biological_unit", ">=", 20, "Sensitivity thresholds: 10 and 50"),
    ]
    for family, metric, op, value, note in rules:
        rows.append({"entry_type": "rule", "marker_family": family, "panel": metric,
                     "gene": None, "operator": op, "threshold": value, "note": note})
    return pd.DataFrame(rows)


def prepare_workbook_data() -> Path:
    BUILD.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    s1 = union_frames([
        ("cohort_overview", read("results/frozen_v1/cohort_overview.csv")),
        ("sample_inclusion", read("results/frozen_v1/sample_inclusion_ledger.csv").drop(columns=["library"], errors="ignore")),
        ("sample_qc", read("results/preflight/sample_qc.csv").drop(columns=["library"], errors="ignore")),
    ])
    s2 = union_frames([
        ("marker_or_rule", marker_table()),
        ("label_transfer_confusion", read("results/advanced_singlecell/cross_dataset_label_transfer_confusion.csv")),
        ("label_transfer_metrics", read("results/advanced_singlecell/cross_dataset_label_transfer_metrics.csv")),
    ])
    s3 = read("results/frozen_v1/formal_composition_effects.csv")

    formal = read("results/frozen_v1/formal_serous_gene_models.csv.gz")
    locked = read("deliverables/stage_03_full_figure_set/tables/locked_gene_model_estimates.csv")[["gene", "locked"]]
    formal = formal.merge(locked, on="gene", how="left")
    formal["locked"] = formal["locked"].fillna(False).astype(bool)
    formal.insert(0, "record_type", "formal_serous_gene_model")
    s4 = formal

    s5 = union_frames([
        ("functional_enrichment", read("results/formal_models/serous_concordant_functional_enrichment.csv")),
        ("regulatory_activity", read("results/advanced_regulatory/regulatory_cross_cohort_summary.csv")),
        ("collectri_target_coverage", read("results/advanced_regulatory/collectri_target_coverage.csv")),
        ("progeny_target_coverage", read("results/advanced_regulatory/progeny_target_coverage.csv")),
    ])
    s6 = union_frames([
        ("bulk_primary_validation", read("results/bulk_crs_validation/serous_module_validation.csv")),
        ("module_size_sensitivity", read("results/bulk_crs_validation/serous_module_sensitivity.csv")),
        ("cell_count_threshold_sensitivity", read("results/anatomy_aware/cell_threshold_sensitivity.csv")),
        ("GSE136825_pair_candidate_audit", gse136825_pairing_audit()),
        ("GSE136825_included_pair_scores", read("results/bulk_crs_validation/GSE136825_paired_module_scores.csv")),
        ("GSE136825_pairing_statistics", read("results/bulk_crs_validation/GSE136825_paired_module_statistics.csv")),
        ("leave_one_donor_out_summary", read("results/robustness/serous_patient_leave_one_out_summary.csv")),
        ("leave_one_donor_out_gene", read("results/robustness/serous_patient_leave_one_out_genes.csv")),
    ])
    s7 = union_frames([
        ("network_summary", read("results/advanced_coexpression/module_preservation_summary.csv")),
        ("network_node", read("results/advanced_coexpression/preserved_module_nodes.csv")),
        ("consensus_edge", read("results/advanced_coexpression/preserved_module_edges.csv")),
        ("unstratified_random_density_500", read("results/advanced_coexpression/module_density_permutation.csv")),
        ("matched_random_density_5000", read("results/advanced_coexpression/module_density_matched_permutation.csv")),
        ("matching_stratum_audit", read("results/advanced_coexpression/module_density_matched_strata_audit.csv")),
    ])
    s8_parts = [
        ("ecm_residualized_score", read("deliverables/section_3_6_boundary_evidence/ecm_coupling_residualized_sample_scores.csv")),
        ("spatial_gene_coverage", read("results/spatial_validation/module_gene_overlap.csv")),
        ("spatial_patient_score", read("results/spatial_validation/patient_module_scores.csv")),
        ("spatial_summary", read("results/spatial_validation/spatial_module_summary.csv")),
        ("protein_coverage", read("results/protein_validation/PXD013330_module_proteins.csv")),
        ("protein_sample_score", read("results/protein_validation/PXD013330_sample_module_scores.csv")),
        ("protein_summary", read("results/protein_validation/PXD013330_module_summary.csv")),
    ]
    s8 = union_frames(s8_parts)

    readme = pd.DataFrame([
        ["Table S1", "Datasets, sample groups, inclusion/exclusion and sample-level QC", "One row per record; use record_type to distinguish sections"],
        ["Table S2", "Predefined marker panels, annotation thresholds, sampling rules and label-transfer audit", "Marker genes were fixed before disease-effect testing"],
        ["Table S3", "All secretory-cell composition comparisons", "Biological unit is donor or individual biopsy, as metadata permit"],
        ["Table S4", "Complete two-cohort serous pseudobulk gene models and locked 40-gene flag", "Definition genes were excluded before gene-effect testing"],
        ["Table S5", "Functional enrichment and CollecTRI/PROGENy activity results", "Transcriptome-derived regulatory scores are mechanistic hypotheses"],
        ["Table S6", "Bulk validation, GSE136825 pairing audit, module-size sensitivity and leave-one-donor-out analyses", "Locked module was not refitted in validation cohorts"],
        ["Table S7", "198-gene co-expression structure, 7,389 consensus edges and random-gene-set nulls", "198 = 38 detectable locked genes + 160 discovery neighbors"],
        ["Table S8", "ECM, GeoMx and secretion-proteomics boundary evidence", "Coverage/direction audit; not clinical biomarker validation"],
    ], columns=["sheet", "contents", "interpretation_note"])

    payload = {
        "title": "Supplementary Data: anatomy-aware multi-cohort analysis of serous glandular dysfunction in CRSwNP",
        "generated": "2026-09-01",
        "sheets": {
            "README": frame_payload(readme, "Workbook guide and interpretation boundaries"),
            "Table S1": frame_payload(s1, "Datasets, sample groups, inclusion/exclusion and QC"),
            "Table S2": frame_payload(s2, "Marker panels, annotation thresholds, integration sampling rules and label-transfer audit"),
            "Table S3": frame_payload(s3, "Secretory-cell composition comparisons"),
            "Table S4": frame_payload(s4, "Complete serous pseudobulk gene models and locked-module membership"),
            "Table S5": frame_payload(s5, "Functional enrichment and inferred regulatory activity"),
            "Table S6": frame_payload(s6, "Bulk validation, GSE136825 pair reconstruction and robustness analyses"),
            "Table S7": frame_payload(s7, "Extended network and random-network analyses"),
            "Table S8": frame_payload(s8, "Orthogonal public-resource boundary evidence"),
        },
    }
    target = BUILD / "supplementary_workbook_data.json"
    target.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    return target


def setup_style():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.8,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })


def panel_label(ax, label):
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


def save_figure(fig, stem):
    fig.savefig(FIG / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def make_s1():
    qc = read("results/preflight/cohort_qc_summary.csv")
    ledger = read("results/frozen_v1/sample_inclusion_ledger.csv")
    audit = read("results/advanced_singlecell/integration_sample_audit.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1.0, 1.35]})
    ax = axes[0]
    stages = ["Readable\nmatrices", "QC-passing\ncells", "Integrated\nsecretory cells"]
    values = [int(ledger.included_readable.sum()), int(qc.cells_ge_200_genes.sum()), int(audit.cells_retained_balanced.sum())]
    y = np.arange(3)
    ax.barh(y, values, color=[NAVY, BLUE, TEAL], height=.55)
    ax.set_yticks(y, stages)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Analysis units (log scale)")
    for yi, v in zip(y, values):
        ax.text(v * 1.08, yi, f"{v:,}", va="center", fontsize=8, color=NAVY)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#DDE3E8", lw=.6)
    ax.set_axisbelow(True)
    panel_label(ax, "A")

    ax = axes[1]
    x = np.arange(len(qc))
    colors = [BLUE if d == "GSE235711" else ORANGE if d == "GSE261706" else TEAL for d in qc.dataset]
    labels = [f"{dis} | {t.replace('_',' ')}" for dis, t in zip(qc.disease, qc.tissue)]
    ax.barh(x, qc.cells_ge_200_genes, color=colors, height=.7)
    ax.set_yticks(x, labels)
    ax.invert_yaxis()
    ax.set_xlabel("QC-passing cells")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#E5E7EB", lw=.6)
    ax.set_axisbelow(True)
    handles = [Line2D([0], [0], color=c, lw=6, label=d) for d, c in [("GSE235711", BLUE), ("GSE261706", ORANGE), ("GSE276503", TEAL)]]
    ax.legend(handles=handles, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(.5, 1.15), fontsize=7)
    panel_label(ax, "B")
    fig.suptitle("Sample and cell-quality audit", x=.5, y=1.03, fontsize=10, fontweight="bold", color=NAVY)
    fig.tight_layout()
    save_figure(fig, "Figure_S1_QC")


def make_s2():
    conf = read("results/advanced_singlecell/cross_dataset_label_transfer_confusion.csv")
    comp = read("results/frozen_v1/formal_composition_effects.csv")
    labels = list(SECRETORY) + ["Ambiguous_secretory"]
    pretty = [x.replace("_", " ") for x in labels]
    fig, axes = plt.subplots(1, 3, figsize=(8.2, 3.2), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    cmap = LinearSegmentedColormap.from_list("journalblue", ["#F7FAFC", "#9EC5E5", "#0E4A7A"])
    for ax, (train, test), lab in zip(axes[:2], [("GSE235711", "GSE276503"), ("GSE276503", "GSE235711")], ["A", "B"]):
        z = conf[(conf.train_dataset == train) & (conf.test_dataset == test)]
        mat = z.pivot(index="truth", columns="predicted", values="fraction").reindex(index=labels, columns=labels).fillna(0)
        im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect="equal")
        ax.set_xticks(range(len(labels)), pretty, rotation=55, ha="right")
        ax.set_yticks(range(len(labels)), pretty)
        ax.set_title(f"{train} → {test}", fontweight="bold", color=NAVY)
        ax.set_xlabel("Transferred identity")
        ax.set_ylabel("Original identity" if lab == "A" else "")
        if lab == "B":
            ax.set_yticklabels([])
        for i in range(len(labels)):
            v = mat.iloc[i, i]
            ax.text(i, i, f"{v:.0%}", ha="center", va="center", color="white" if v > .45 else NAVY, fontsize=6.5, fontweight="bold")
        for s in ax.spines.values(): s.set_visible(False)
        panel_label(ax, lab)

    ax = axes[2]
    z = comp[(comp.endpoint_role == "secondary_subtype") & (comp.dataset == "GSE276503") &
             (comp.contrast == "NP_vs_CRSwNP_IT")].copy()
    z = z.sort_values("median_difference")
    yy = np.arange(len(z))
    ax.hlines(yy, z.ci95_low, z.ci95_high, color="#9AA5B1", lw=1.3)
    ax.scatter(z.median_difference, yy, c=np.where(z.median_difference < 0, BLUE, ORANGE), s=32, zorder=3, edgecolor="white", linewidth=.5)
    ax.axvline(0, color="#4B5563", lw=.8, ls="--")
    ax.set_yticks(yy, [x.replace("_", " ") for x in z.subtype])
    ax.set_xlabel("Median fraction difference\n(CRSwNP polyp − non-polyp mucosa)")
    ax.set_title("Secondary composition effects", fontweight="bold", color=NAVY)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#E5E7EB", lw=.6)
    panel_label(ax, "C")
    fig.suptitle("Cross-dataset annotation transfer and secondary composition endpoints", y=1.02, fontsize=10, fontweight="bold", color=NAVY)
    fig.subplots_adjust(left=.09, right=.98, bottom=.28, top=.82, wspace=.45)
    save_figure(fig, "Figure_S2_LabelTransfer_Composition")


def make_s3():
    sens = read("results/anatomy_aware/cell_threshold_sensitivity.csv")
    loo = read("results/robustness/serous_patient_leave_one_out_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ax = axes[0]
    types = list(dict.fromkeys(sens.celltype))
    markers = {20: "o", 50: "s"}
    type_colors = dict(zip(types, [NAVY, TEAL, ORANGE, RED]))
    for ct, color in type_colors.items():
        z = sens[sens.celltype == ct]
        for row in z.itertuples():
            ax.scatter(row.rho, row.same_direction_rate, marker=markers.get(row.cell_min, "o"), s=48, color=color, edgecolor="white", linewidth=.6)
    ax.set_xlabel("Cross-cohort Spearman ρ")
    ax.set_ylabel("Directional agreement")
    ax.set_ylim(.48, .84)
    ax.set_xlim(-.02, .36)
    ax.yaxis.set_major_formatter(lambda x, pos: f"{x:.0%}")
    ax.grid(color="#E5E7EB", lw=.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=c, label=ct)
                       for ct, c in type_colors.items()], frameon=False, loc="upper left", ncol=2, fontsize=6.5,
              handletextpad=.3, columnspacing=.8)
    ax.text(.98, .04, "circle: ≥20 cells   square: ≥50 cells", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.3, color=GREY)
    ax.set_title("Cell-count threshold sensitivity", fontweight="bold", color=NAVY)
    panel_label(ax, "A")

    ax = axes[1]
    z = loo.copy()
    labels = z.omitted_donor.replace({"none_full": "Full"}).str.replace("CRSwNP_", "−P", regex=False)
    y = np.arange(len(z))
    ax.hlines(y, 0, z.rho_external, color="#B7C2CC", lw=2)
    ax.scatter(z.rho_external, y, s=45, color=np.where(z.omitted_donor == "none_full", ORANGE, BLUE), edgecolor="white", linewidth=.7, zorder=3)
    for yi, rr in zip(y, z.rho_external):
        ax.text(rr + .008, yi, f"{rr:.3f}", va="center", fontsize=7, color=NAVY)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, .36)
    ax.set_xlabel("Correlation with external-cohort effects")
    ax.set_title("Leave-one-paired-donor-out stability", fontweight="bold", color=NAVY)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#E5E7EB", lw=.6)
    panel_label(ax, "B")
    fig.suptitle("Robustness to cell-count thresholds and donor influence", y=1.04, fontsize=10, fontweight="bold", color=NAVY)
    fig.tight_layout()
    save_figure(fig, "Figure_S3_Robustness")


def make_s4():
    audit = read("results/advanced_coexpression/module_density_matched_strata_audit.csv")
    null = read("results/advanced_coexpression/module_density_matched_permutation.csv").iloc[:, 0]
    summary = read("results/advanced_coexpression/module_preservation_summary.csv").iloc[0]
    grid = np.zeros((5, 5))
    for r in audit.itertuples():
        i, j = map(int, r.stratum.split("_"))
        grid[j, i] = r.module_genes
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1, 1.35]})
    ax = axes[0]
    im = ax.imshow(grid, origin="lower", cmap=LinearSegmentedColormap.from_list("teal", ["#F5FAFA", "#66B5AE", "#087F74"]), vmin=0)
    for j in range(5):
        for i in range(5):
            if grid[j, i] > 0:
                ax.text(i, j, str(int(grid[j, i])), ha="center", va="center", fontsize=8, color="white" if grid[j, i] > 20 else NAVY, fontweight="bold")
    ax.set_xticks(range(5), ["Q1", "Q2", "Q3", "Q4", "Q5"])
    ax.set_yticks(range(5), ["Q1", "Q2", "Q3", "Q4", "Q5"])
    ax.set_xlabel("Mean-expression quintile")
    ax.set_ylabel("Residual-variance quintile")
    ax.set_title("Exact matched-stratum counts", fontweight="bold", color=NAVY)
    for s in ax.spines.values(): s.set_visible(False)
    panel_label(ax, "A")

    ax = axes[1]
    ax.hist(null, bins=42, color="#A8C8DE", edgecolor="white", linewidth=.35, density=True)
    observed = float(summary.validation_absolute_density)
    ax.axvline(observed, color=RED, lw=2)
    ax.text(observed, ax.get_ylim()[1] * .82, f"Observed = {observed:.3f}\nEmpirical p = {summary.validation_density_matched_empirical_p:.4g}",
            ha="right", va="top", fontsize=8, color=RED)
    ax.set_xlabel("Validation absolute co-expression density")
    ax.set_ylabel("Density")
    ax.set_title("5,000 expression/variance-matched gene sets", fontweight="bold", color=NAVY)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E5E7EB", lw=.6)
    panel_label(ax, "B")
    fig.suptitle("Matched random-network sensitivity analysis", y=1.04, fontsize=10, fontweight="bold", color=NAVY)
    fig.tight_layout()
    save_figure(fig, "Figure_S4_MatchedNetwork")


def make_s5():
    ecm = read("deliverables/section_3_6_boundary_evidence/ecm_coupling_residualized_sample_scores.csv")
    coverage = read("results/spatial_validation/module_gene_overlap.csv")
    spatial = read("results/spatial_validation/patient_module_scores.csv")
    protein = read("results/protein_validation/PXD013330_sample_module_scores.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.7))
    ax = axes[0, 0]
    for ds, color, marker in [("GSE36830", BLUE, "o"), ("GSE136825", TEAL, "s")]:
        z = ecm[ecm.dataset == ds]
        ax.scatter(z.serous_residual, z.ecm_residual, s=26, color=color, alpha=.75, edgecolor="white", linewidth=.4, marker=marker, label=ds)
        if len(z) > 1:
            co = np.polyfit(z.serous_residual, z.ecm_residual, 1)
            xx = np.linspace(z.serous_residual.min(), z.serous_residual.max(), 50)
            ax.plot(xx, co[0] * xx + co[1], color=color, lw=1.2)
    ax.axhline(0, color="#CDD5DD", lw=.7); ax.axvline(0, color="#CDD5DD", lw=.7)
    ax.set_xlabel("Residualized serous-module score")
    ax.set_ylabel("Residualized ECM score")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("ECM coupling", fontweight="bold", color=NAVY)
    ax.spines[["top", "right"]].set_visible(False)
    panel_label(ax, "A")

    ax = axes[0, 1]
    coverage = coverage.sort_values(["available", "gene"], ascending=[False, True]).reset_index(drop=True)
    x = np.arange(len(coverage))
    ax.scatter(x, np.zeros(len(x)), s=np.where(coverage.available, 35, 12), c=np.where(coverage.available, TEAL, "#D6DEE5"), marker="s")
    ax.set_xlim(-1, 40); ax.set_ylim(-.25, .65)
    ax.set_yticks([]); ax.set_xticks([])
    ax.text(.98, .86, f"{int(coverage.available.sum())}/40 genes available", transform=ax.transAxes, ha="right", fontsize=9, fontweight="bold", color=TEAL)
    available_names = ", ".join(coverage.loc[coverage.available, "gene"].astype(str))
    ax.text(.02, .64, f"Available genes: {available_names}", transform=ax.transAxes,
            ha="left", va="center", fontsize=7.5, color=NAVY)
    ax.set_title("GeoMx panel coverage", fontweight="bold", color=NAVY)
    for s in ax.spines.values(): s.set_visible(False)
    panel_label(ax, "B")

    ax = axes[1, 0]
    rng = np.random.default_rng(SEED)
    groups = [("PanCK", "CTRL"), ("PanCK", "CRSwNP"), ("CD45", "CTRL"), ("CD45", "CRSwNP")]
    xpos = [0, 1, 2.6, 3.6]
    for x0, (comp, grp) in zip(xpos, groups):
        z = spatial[(spatial.compartment == comp) & (spatial.group == grp)].score.to_numpy()
        jit = rng.uniform(-.08, .08, len(z))
        color = NAVY if grp == "CTRL" else ORANGE
        ax.scatter(np.full(len(z), x0) + jit, z, s=34, color=color, edgecolor="white", linewidth=.6, zorder=3)
        if len(z): ax.hlines(np.median(z), x0 - .22, x0 + .22, color="#111827", lw=1.5)
    ax.set_xticks(xpos, ["Healthy", "CRSwNP", "Healthy", "CRSwNP"], rotation=25, ha="right")
    ax.text(.5, -.25, "PanCK", transform=ax.get_xaxis_transform(), ha="center", fontsize=8, fontweight="bold")
    ax.text(3.1, -.25, "CD45", transform=ax.get_xaxis_transform(), ha="center", fontsize=8, fontweight="bold")
    ax.set_ylabel("Four-gene spatial proxy")
    ax.set_title("Patient-level spatial proxy", fontweight="bold", color=NAVY)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", color="#E5E7EB", lw=.6)
    panel_label(ax, "C")

    ax = axes[1, 1]
    order = ["CON", "CRSwNP"]
    for x0, grp in enumerate(order):
        z = protein[protein.group == grp].score.to_numpy()
        jit = rng.uniform(-.08, .08, len(z))
        color = NAVY if grp == "CON" else ORANGE
        ax.scatter(np.full(len(z), x0) + jit, z, s=42, color=color, edgecolor="white", linewidth=.7)
        ax.hlines(np.median(z), x0 - .22, x0 + .22, color="#111827", lw=1.5)
    ax.set_xticks([0, 1], ["Control", "CRSwNP"])
    ax.set_ylabel("Seven-protein proxy score")
    ax.text(.5, .03, "Pooled biological material; technical replicates", transform=ax.transAxes, ha="center", fontsize=6.5, color=GREY)
    ax.set_title("Nasal-secretion proteomics", fontweight="bold", color=NAVY)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", color="#E5E7EB", lw=.6)
    panel_label(ax, "D")
    fig.suptitle("Orthogonal public resources define translation boundaries", y=.99, fontsize=10, fontweight="bold", color=NAVY)
    fig.subplots_adjust(left=.1, right=.97, bottom=.12, top=.9, hspace=.52, wspace=.35)
    save_figure(fig, "Figure_S5_BoundaryEvidence")


def make_figures():
    setup_style()
    FIG.mkdir(parents=True, exist_ok=True)
    make_s1(); make_s2(); make_s3(); make_s4(); make_s5()


def build_pdf():
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER, TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, PageBreak,
                                    Image, KeepTogether, Table, TableStyle)

    target = OUT / "Supplementary Methods and Figures.pdf"
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="SupTitle", parent=styles["Title"], fontName="Helvetica-Bold", fontSize=17,
                              leading=21, textColor=colors.HexColor(NAVY), alignment=TA_CENTER, spaceAfter=12))
    styles.add(ParagraphStyle(name="SupH1", parent=styles["Heading1"], fontName="Helvetica-Bold", fontSize=12,
                              leading=15, textColor=colors.HexColor(NAVY), spaceBefore=10, spaceAfter=6))
    styles.add(ParagraphStyle(name="SupH2", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=9.5,
                              leading=12, textColor=colors.HexColor(BLUE), spaceBefore=7, spaceAfter=3))
    styles.add(ParagraphStyle(name="SupBody", parent=styles["BodyText"], fontName="Helvetica", fontSize=8.4,
                              leading=12, alignment=TA_LEFT, spaceAfter=6))
    styles.add(ParagraphStyle(name="Legend", parent=styles["BodyText"], fontName="Helvetica", fontSize=7.8,
                              leading=10.8, spaceAfter=8))
    styles.add(ParagraphStyle(name="Small", parent=styles["BodyText"], fontName="Helvetica", fontSize=7,
                              leading=9, textColor=colors.HexColor(GREY)))

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor("#D7DEE5")); canvas.line(20*mm, 13*mm, 190*mm, 13*mm)
        canvas.setFont("Helvetica", 7); canvas.setFillColor(colors.HexColor(GREY))
        canvas.drawString(20*mm, 8*mm, "Supplementary Methods and Figures")
        canvas.drawRightString(190*mm, 8*mm, str(doc.page))
        canvas.restoreState()

    doc = SimpleDocTemplate(str(target), pagesize=A4, leftMargin=20*mm, rightMargin=20*mm,
                            topMargin=18*mm, bottomMargin=18*mm, title="Supplementary Methods and Figures")
    story = []
    title = "Anatomy-aware multi-cohort transcriptomics distinguishes serous glandular cell depletion from cell-intrinsic host-defense suppression in chronic rhinosinusitis with nasal polyps"
    story += [Spacer(1, 12*mm), Paragraph("Supplementary Methods and Figures", styles["SupTitle"]),
              Paragraph(title, ParagraphStyle(name="ArticleTitle", parent=styles["SupBody"], fontSize=11, leading=15,
                                               alignment=TA_CENTER, textColor=colors.HexColor(NAVY))),
              Spacer(1, 6*mm), Paragraph("Jingren Yan and Shun Ding", ParagraphStyle(name="Authors", parent=styles["SupBody"], alignment=TA_CENTER)),
              Spacer(1, 8*mm)]
    box = Table([[Paragraph("Contents", styles["SupH2"]), Paragraph("Supplementary Methods; Supplementary Figures S1–S5; Supplementary Figure Legends; Supplementary References", styles["SupBody"])],
                 [Paragraph("Companion file", styles["SupH2"]), Paragraph("Supplementary Data.xlsx (Tables S1–S8)", styles["SupBody"])]], colWidths=[38*mm, 125*mm])
    box.setStyle(TableStyle([("BACKGROUND", (0,0), (-1,-1), colors.HexColor("#F5F8FA")), ("BOX", (0,0), (-1,-1), .6, colors.HexColor("#CBD5DF")),
                                  ("INNERGRID", (0,0), (-1,-1), .3, colors.HexColor("#E1E7EC")), ("VALIGN", (0,0), (-1,-1), "TOP"),
                                  ("LEFTPADDING", (0,0), (-1,-1), 7), ("RIGHTPADDING", (0,0), (-1,-1), 7)]))
    story += [box, PageBreak(), Paragraph("Supplementary Methods", styles["SupH1"])]

    methods = [
        ("1. Study design and evidence hierarchy",
         "This secondary analysis used two single-cell RNA-sequencing discovery/replication cohorts (GSE235711 and GSE276503), two independent bulk-transcriptomic validation cohorts (GSE36830 and GSE136825), and public GeoMx and nasal-secretion proteomic resources as boundary evidence. Donor–tissue combinations were the biological unit when donor identifiers were available; otherwise, individual biopsies were treated as independent units. Cells were never used as inferential replicates. Cell-composition analyses, within-serous-cell pseudobulk analyses, locked-module validation, regulatory scoring, cross-cohort co-expression reproducibility, and orthogonal boundary audits were kept analytically distinct."),
        ("2. Cell-level quality control and broad annotation",
         "Raw count matrices were oriented against their feature tables and collapsed to gene symbols when required. Cells were retained when at least 200 genes were detected, total library size was greater than zero, and the mitochondrial fraction was below 20%. Counts used for annotation were normalized per cell to 10,000 and transformed as log(1+x). Broad compartment scores were the mean normalized expression of predefined epithelial, immune, fibroblast, and endothelial marker panels (Table S2). A cell was eligible for epithelial subtyping only when the epithelial panel was the highest-scoring broad compartment, at least one epithelial marker was detected, and the epithelial score was at least 0.15."),
        ("3. Secretory-cell panels and deterministic subtype assignment",
         "Disease-independent marker panels were fixed before disease-effect testing: Surface/club (SCGB1A1, SCGB3A1, KRT4, KRT13, KLF3, CYP2F1); Serous glandular (LTF, LYZ, BPIFB1, SLPI, WFDC2, PRR4, LPO); Mucous glandular (MUC5B, BPIFB2, AGR2, SPDEF, ZG16B, PIP); Goblet (MUC5AC, CLCA1, SPDEF, AGR2, FCGBP); Inflammatory secretory (CST1, SLC26A4, SERPINB3, SERPINB4, CCL26, POSTN); and Transitional secretory (KRT4, KRT13, KRT17, KRT19, SFN, KRT16). The winning subtype required at least one detected marker and a mean score of at least 0.12. Cells with a best-minus-second-best score margin below 0.05 were labeled Ambiguous_secretory. Immunoglobulin fractions were retained as diagnostics and were not automatic exclusion criteria."),
        ("4. Balanced cross-cohort integration of 36,116 cells",
         "To avoid domination by large samples or abundant states, at most 350 cells were sampled without replacement from each sample-by-subtype stratum using random seed 20260831. This deterministic balanced sample contained 36,116 secretory cells. The common-gene feature space excluded all annotation-definition genes, mitochondrial genes, ribosomal genes (RPL/RPS), immunoglobulin genes, and genes with insufficient variance. The 2,500 highest-variance eligible genes were used, supplemented only by required audit features when shared between cohorts. Truncated singular-value decomposition generated 40 components, followed by Harmony correction for dataset. The first 30 corrected components supported visualization and reciprocal distance-weighted 25-nearest-neighbor label transfer. t-SNE (perplexity 50, Barnes–Hut, 1,200 iterations) and eight-cluster k-means were descriptive only and were not used for statistical inference."),
        ("5. Pseudobulk construction, gene filtering and 40-gene module locking",
         "Raw counts were summed within each donor–tissue–secretory-subtype unit and converted to log2(CPM+0.5). The primary eligibility threshold was at least 20 cells per biological unit, with 10- and 50-cell thresholds used for sensitivity analyses. Genes used to define broad or secretory identities were removed before testing cell-intrinsic effects, as were genes detected in fewer than two combined biological units. GSE235711 used paired nasal-polyp minus ethmoid effects; GSE276503 used unpaired nasal-polyp minus inferior-turbinate effects. Candidate module genes had effects of at most −1 log2 units in both cohorts. Vascular, stromal, erythroid, immunoglobulin, HLA, ribosomal, and non-canonical feature symbols were excluded. Remaining genes were ranked by the smaller absolute cross-cohort effect and the leading 40 were locked before bulk validation. Formal gene-wise estimates included cohort-specific standard errors, confidence intervals, P values, Benjamini–Hochberg FDR, fixed-effect synthesis, Cochran Q and I²."),
        ("6. Bulk-cohort preprocessing and GSE136825 pair reconstruction",
         "For GSE136825, the featureCounts sample columns were selected after the five annotation columns. Ensembl version suffixes were removed, identifiers were mapped to gene symbols using the downloaded GSE235711 10x feature table, and rows mapping to the same symbol were summed. Library-size-normalized counts per million were calculated and transformed as log2(CPM+0.5) in the analysis code. Module scores were then computed as the mean of gene-wise standardized expression values across the detected locked genes. Clinical groups were assigned from the GEO !Sample_source_name_ch1 field. The paired analysis did not use BAM sample names or a separate patient metadata field: pairs were reconstructed from GEO !Sample_title. Titles PY <identifier> and SCP <identifier> formed one series, and TT <identifier> and SCT <identifier> formed the second series; only records with the same series and identical numeric suffix were merged. All role/series/identifier keys were unique, so no duplicate-resolution procedure was invoked. The PY/SCP series contained 29 polyp and 19 inferior-turbinate candidates, yielding 17 matched pairs; the TT/SCT series contained 13 and 14 candidates, yielding 13 matched pairs. Thus, 30 pairs were retained, while 12 polyp and three inferior-turbinate candidates without a valid counterpart were excluded. The complete candidate-level audit and the 30 retained pairs are provided in Table S6. The main-text wording 'log2-transformed counts per million' can remain unchanged; this implementation detail does not require an additional citation."),
        ("7. Serous-state diffusion coordinate",
         "A serous-only geometry was recomputed using the 1,500 most variable eligible non-definition genes. Thirty truncated-SVD components were followed by Harmony correction for dataset. A 25-nearest-neighbor Euclidean graph was formed on the first 25 corrected components, with affinities exp(−d²/σ²), where σ was the median 25th-neighbor distance. The symmetrized affinity matrix was degree-normalized and decomposed with a seeded ARPACK start vector. The trivial eigenvector was removed and ten diffusion coordinates were retained. The root was the healthy GSE235711 serous cell maximizing locked-module score minus cilia score. Euclidean distance from this root in diffusion space was converted to a percentile-ranked serous-state coordinate. Donor/tissue summaries used the median cell value. This coordinate is descriptive and is not interpreted as pseudotime, lineage direction, or functional failure."),
        ("8. Functional enrichment and regulatory activity inference",
         "Functional enrichment was performed separately for concordantly decreased and increased non-definition genes; full tested terms, overlap genes, nominal P values and within-library FDR are supplied in Table S5. CollecTRI and PROGENy resources were frozen before cohort contrasts. For PROGENy, annotations with P<0.05 were retained and the 500 strongest responsive genes per pathway were used. CollecTRI regulons required at least 20 measured targets and PROGENy pathways at least 25. Within each cohort, gene-wise logCPM values were standardized and activity was calculated as the signed weighted sum divided by the L2 norm of network weights. GSE235711 used paired one-sample contrasts with leave-one-pair-out sign stability; GSE276503 used Welch contrasts. Inverse-variance fixed-effect estimates, Fisher-combined P values and Benjamini–Hochberg FDR were reported. These scores are transcriptome-derived hypotheses and do not demonstrate causal regulator activity."),
        ("9. Construction and cross-cohort reproducibility of the 198-gene co-expression structure",
         "Bulk expression was residualized within dataset-specific clinical groups. Of the locked 40 genes, 38 were detectable in the common GSE136825/GSE36830 expression space. A standardized 38-gene score was computed in the discovery cohort GSE136825. After excluding locked genes, the 160 genes with the highest positive discovery Pearson correlation with this score were frozen, producing a 198-gene module. All 19,503 gene pairs were used for cross-cohort comparison of within-cohort Pearson coefficients and absolute-density statistics; their cross-cohort rank association was summarized by Spearman correlation. A readable consensus network retained positive edges with discovery Pearson correlation at least 0.55 and validation Pearson correlation at least 0.30; consensus weight was the geometric mean of the two correlations. This yielded 7,389 edges. GSE36830 served as the validation cohort and was not used to select neighbors or thresholds."),
        ("10. Random-network analyses",
         "The initial null comprised 500 unstratified same-size gene sets drawn without replacement. The stronger prespecified sensitivity null comprised 5,000 same-size sets exactly matched to the 198-gene module across a 5×5 rank-quantile grid of validation-cohort mean expression and within-group residual variance. Module genes were excluded from the eligible non-module universe (n=18,052), and sampling occurred without replacement within each stratum. Each null set was scored as the mean absolute validation-cohort pairwise correlation. Empirical P values used (1 + number of null statistics at least as large as observed)/(1 + number of permutations). Random seed was 20260831."),
        ("11. Orthogonal public-resource boundary analyses",
         "Within GSE36830 and GSE136825, serous-module and 11-gene extracellular-matrix scores were residualized for dataset-specific group effects before Spearman correlation. GeoMx coverage was limited to four of 40 locked genes; patient-level four-gene proxy scores were therefore treated as coverage and direction audits in PanCK- and CD45-enriched regions. PXD013330 contained eleven detected module-mapped proteins, of which seven were quantified across all nine MaxQuant columns. Because these columns were technical replicates of pooled material and the specimen was nasal secretion rather than tissue, the seven-protein score was not interpreted as patient-level validation or a non-invasive biomarker test."),
        ("12. Software, reproducibility and code availability",
         "Core analyses were run with Python 3.12.13, NumPy 2.4.6, pandas 3.0.5, SciPy 1.15.3, scikit-learn 1.7.1, anndata 0.12.2, harmonypy 0.0.10 and decoupler 2.1.1. Publication figures were generated with R 4.6.1, ggplot2 4.0.3, ragg 1.5.2 and scales 1.4.0; selected supplementary composites were assembled with Matplotlib 3.11.1. The global stochastic seed was 20260831. Intermediate tables were written after each analytical stage, and all reported values are traceable to Supplementary Data.xlsx. Analysis code is available from the authors during review and should be deposited in a versioned public repository before publication; the repository DOI/URL should replace this statement at submission."),
    ]
    for heading, body in methods:
        story.append(KeepTogether([Paragraph(heading, styles["SupH2"]), Paragraph(body, styles["SupBody"])]))

    story += [PageBreak(), Paragraph("Supplementary Figures", styles["SupH1"])]
    legends = {
        "Figure S1": ("Figure_S1_QC.png", "Sample and cell-quality audit. (A) Readable matrices, QC-passing cells and the deterministically balanced 36,116-cell integration set. The logarithmic axis is used only to display different count scales. (B) QC-passing cells by dataset, disease and anatomical site. Bars summarize cells meeting the predefined detected-gene threshold; sample-level mitochondrial and library metrics are provided in Table S1."),
        "Figure S2": ("Figure_S2_LabelTransfer_Composition.png", "Cross-dataset annotation transfer and secondary composition endpoints. (A,B) Complete reciprocal distance-weighted 25-nearest-neighbor transfer matrices in Harmony-corrected space. Percentages are shown on the diagonal; the full matrices and global balanced accuracy/macro-F1 are supplied in Supplementary Data. (C) Secondary secretory-subtype composition effects shown as median fraction differences with 95% confidence intervals. These endpoints support annotation and composition auditing and are not treated as independent cell-level tests."),
        "Figure S3": ("Figure_S3_Robustness.png", "Robustness to cell-count thresholds and donor influence. (A) Cross-cohort correlation and directional agreement across available cell types under minimum 20- and 50-cell pseudobulk thresholds. (B) Correlation between paired-cohort and external-cohort serous gene effects in the full analysis and after omitting each paired donor. All 40 locked genes retained negative paired effects in every leave-one-donor-out iteration (Table S6)."),
        "Figure S4": ("Figure_S4_MatchedNetwork.png", "Expression- and variance-matched random-network analysis. (A) Exact module-gene counts across the 5×5 grid of validation-cohort mean-expression and residual-variance quintiles. Empty strata contained no module genes and therefore required no matched draws. (B) Validation absolute co-expression density for 5,000 matched random 198-gene sets. The red line marks the observed validation density (0.379); the empirical P value was 0.00020."),
        "Figure S5": ("Figure_S5_BoundaryEvidence.png", "Orthogonal public resources define translation boundaries. (A) Residualized serous-module and 11-gene ECM scores in two bulk cohorts. Neither cohort showed significant coupling. (B) Only four of 40 locked genes were represented on the available GeoMx panel. (C) Patient-level four-gene spatial proxy scores in PanCK- and CD45-enriched compartments; horizontal lines denote medians. (D) Seven-protein nasal-secretion proxy across MaxQuant technical-replicate columns from pooled biological material. These analyses are coverage/direction audits and do not constitute spatial or non-invasive biomarker validation."),
    }
    for label, (fname, legend) in legends.items():
        img = Image(str(FIG / fname), width=168*mm, height=0)
        iw, ih = img.imageWidth, img.imageHeight
        img.drawWidth = 168*mm
        img.drawHeight = img.drawWidth * ih / iw
        story += [KeepTogether([Paragraph(label, styles["SupH2"]), img, Spacer(1, 2*mm), Paragraph(f"<b>{label}.</b> {legend}", styles["Legend"])])]
        if label in {"Figure S2", "Figure S4"}: story.append(PageBreak())

    story += [PageBreak(), Paragraph("Supplementary References", styles["SupH1"])]
    refs = [
        "1. Korsunsky I, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nature Methods. 2019;16:1289–1296.",
        "2. Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances. 2022;2:vbac016.",
        "3. Schubert M, et al. Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications. 2018;9:20.",
        "4. Müller-Dott S, et al. Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities. Nucleic Acids Research. 2023;51:10934–10949.",
    ]
    for r in refs: story.append(Paragraph(r, styles["SupBody"]))
    story += [Spacer(1, 4*mm), Paragraph("Important submission note: replace the provisional code-availability sentence with a permanent repository URL or DOI before journal submission.", styles["Small"])]
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    return target


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["prepare", "pdf"])
    args = ap.parse_args()
    if args.mode == "prepare":
        path = prepare_workbook_data()
        make_figures()
        print(path)
    else:
        print(build_pdf())


if __name__ == "__main__":
    main()
