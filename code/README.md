# Code inventory and execution map

All scripts in this directory were used during the project. They retain their actual filenames; no artificial numbered wrappers were created. Run scripts from the repository root.

## Core single-cell preparation and annotation

1. `prepare_and_qc.py` — archive extraction, public-sample manifest, and fixed QC audit
2. `preflight_signatures.py` — initial marker/signature feasibility checks
3. `celltype_preflight.py` — broad-compartment annotation using predefined panels
4. `secretory_subtype_audit.py` — disease-independent secretory-state assignment
5. `build_pseudobulk.py` — biological-unit-level pseudobulk construction
6. `secretory_subtype_pseudobulk_validation.py` — non-definition-gene subtype audit

## Composition and within-cell transcriptional analysis

- `patient_level_effects.py`
- `pseudobulk_effects.py`
- `anatomy_aware_audit.py`
- `secretory_composition_effects.py`
- `build_frozen_patient_objects.py`
- `formal_serous_pseudobulk_models.py`
- `formal_function_and_site_models.py`

These scripts generate unit-level composition estimates, serous pseudobulks, cross-cohort effects, anatomy-aware contrasts, and functional-enrichment inputs. Cells are not treated as inferential replicates.

## Integration, state coordinate, and regulatory context

- `advanced_singlecell_integration.py` — Harmony integration, t-SNE display, reciprocal 25-NN label transfer, and diffusion-state calculation
- `recompute_serous_continuum_deterministic.py` — deterministic ARPACK `v0` recomputation
- `audit_deterministic_state_recompute.py` — comparison with archived state coordinates
- `regulatory_activity_audit.py` — frozen CollecTRI and PROGENy activity scoring

## External bulk validation and robustness

- `validate_serous_bulk_crs.py` — locked-score assessment in GSE36830 and GSE136825
- `compute_gse136825_paired_module_stats.R` — exact paired GSE136825 module statistics
- `serous_module_sensitivity_coupling.py` — sensitivity and ECM coupling analyses
- `serous_patient_leave_one_out.py` — leave-one-patient-out stability

The locked module is read from `results/locked_40_gene_module.csv`. The public repository uses evidence-bounded terminology consistently; no functional-failure claim is encoded in the filename.

## Cross-cohort co-expression reproducibility

- `bulk_coexpression_preservation.py` — freezes the 198-gene structure, calculates all within-cohort Pearson gene-pair coefficients, consensus edges, and random-gene-set nulls
- `export_all_coexpression_pairs.py` — exports all 19,503 Pearson pairs and their cross-cohort Spearman concordance

The word `preservation` remains in one historical filename for provenance. The manuscript terminology is **cross-cohort co-expression reproducibility**, not formal network preservation.

## Exploratory boundary audits

- `analyze_spatial_validation.py` — four-gene GeoMx proxy analysis
- `reanalyze_pxd013330_proteome.py` — descriptive pooled technical-column proteomic proxy

These analyses are exploratory and do not independently validate the complete module.

## Figure generation

- Figure 1: `plot_article_workflow_clean_final.R`
- Figure 2: `plot_figure2_highimpact_v2.R`
- Figure 3: `export_figure3_final_tiffs.R`, with `plot_figure3_editorial_v4.R`, `plot_figure3_intrinsic_v1.R`, and the targeted `rebuild_fig3b_lowercase_p.R`
- Figure 4: `export_figure4_final_90mm.R`
- Figure 5: `export_figure5_final_90mm.R`
- Figure 6: `export_figure6_final_90mm.R`; matched sensitivity: `export_figure6_matched_sensitivity_90mm.R`
- Figure 7: `export_figure7_highimpact_v2.R`
- Supplementary package: `build_supplementary_package.py`

## Important dependency notes

- `build_pseudobulk.py` imports marker definitions from `celltype_preflight.py`.
- `secretory_subtype_pseudobulk_validation.py` imports panels from `secretory_subtype_audit.py`.
- `formal_serous_pseudobulk_models.py`, `regulatory_activity_audit.py`, and `serous_patient_leave_one_out.py` reuse the frozen pseudobulk loader.
- Co-expression scripts reuse bulk parsing functions from `validate_serous_bulk_crs.py`.
- Figure 3's final exporter sources its two plotting helper scripts.

Exact parameters and package versions are in `environment/analysis_parameters_and_software.md`.
