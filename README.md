# Anatomy-aware multi-cohort transcriptomic analysis of CRSwNP

This repository contains the analysis code and frozen processed outputs accompanying the manuscript:

> *Anatomy-aware multi-cohort transcriptomics distinguishes reduced serous-glandular representation from transcriptional host-defense suppression in chronic rhinosinusitis with nasal polyps.*

The study separates two analytically distinct findings: variation in the relative representation of serous-glandular cells and reproducible within-cell suppression of a serous host-defense transcriptional program. It does not claim absolute cell depletion, a clinically validated biomarker, causal dysfunction, or treatment-response prediction.

## Public datasets

- GSE235711: paired single-cell RNA-sequencing cohort
- GSE276503: external multi-site single-cell RNA-sequencing cohort
- GSE36830: external bulk-transcriptomic validation cohort
- GSE136825: external bulk-transcriptomic validation cohort with a paired polyp–inferior-turbinate subset
- GSE235714: exploratory GeoMx spatial proxy
- PXD013330: exploratory nasal-secretion proteomics

Raw public data are not redistributed. Obtain them from GEO or ProteomeXchange and place them under `data/raw/` as described in [DATA_ACCESS.md](DATA_ACCESS.md).

## Repository structure

- `code/`: scripts actually used in the analysis and final figure generation
- `metadata/`: dataset roles, public sample grouping, pairing audit, and predefined marker panels
- `results/`: frozen processed outputs required to audit or recreate the principal analyses
- `environment/`: historical environment audit and proposed clean reproduction environment
- `deliverables/`: compact plotting inputs used by the final R figure scripts
- `supplementary/`: final Supplementary Data workbook and Supplementary Methods/Figures PDF

Large single-cell objects, raw matrices, local dependency folders, and rendered figure files are deliberately excluded.

## Main analysis workflow

1. Fixed quality control and disease-independent cell annotation
2. Cross-cohort integration and reciprocal label transfer
3. Biological-unit-level composition analysis
4. Serous-glandular pseudobulk analysis with definition genes excluded
5. Construction and locking of the 40-gene serous-defense module
6. Independent bulk-cohort evaluation without refitting
7. Functional enrichment, diffusion-state, CollecTRI, and PROGENy analyses
8. Cross-cohort co-expression reproducibility and random-gene-set sensitivity analyses
9. Exploratory spatial, ECM, and secretion-proteomic boundary audits
10. Main and supplementary figure generation

The exact script inventory and execution dependencies are documented in [code/README.md](code/README.md). Scripts are run from the repository root so their relative paths resolve consistently.

## Key reproducibility details

- Random seed: `20260831` for all explicitly stochastic analyses.
- Biological unit: donor when participant identifiers were available; biopsy otherwise.
- Definition genes used for annotation were excluded before within-serous-cell transcriptional testing.
- The locked 40-gene module was defined using the two single-cell cohorts and assessed in bulk cohorts without refitting.
- The 198-gene co-expression structure contains 38 platform-detectable locked genes plus 160 neighbors selected only in GSE136825.
- Within-cohort gene–gene coefficients are Pearson correlations. Cross-cohort concordance across all 19,503 Pearson coefficients is summarized by Spearman correlation.
- Consensus edges require discovery Pearson `r >= 0.55` and validation Pearson `r >= 0.30`.

Full formulas, thresholds, randomization universes, and software versions are recorded in [environment/analysis_parameters_and_software.md](environment/analysis_parameters_and_software.md).

## Reproducing figures from frozen outputs

Install the recorded R packages, start R in the repository root, and source the relevant figure script. For example:

```r
source("code/export_figure6_final_90mm.R")
```

The scripts write newly rendered files under `deliverables/`. Final publication graphics were exported with Arial, physical dimensions specified in millimetres, and 600-dpi LZW TIFF output; readability was checked at final physical size rather than inferred from dpi alone.

## Environments

For a clean Python reproduction environment:

```bash
conda env create -f environment/environment-unified-python.yml
conda activate crswnp-serous-repro
```

The project used historically distinct early-classification and advanced-analysis environments. They are reported separately and are not retrospectively misrepresented as one exact historical lock. See [environment/README_environment_reproducibility.md](environment/README_environment_reproducibility.md).

## Citation

Please cite the associated manuscript when using this code. Citation details will be updated after publication.

## License

Code is released under the [MIT License](LICENSE). Public-source data and third-party resources remain subject to their original licenses and terms.

## Contact

Jingren Yan  
Department of Medicine 1, Universitätsklinikum Erlangen  
Friedrich-Alexander-Universität Erlangen-Nürnberg  
Email: yanjingren97@163.com
