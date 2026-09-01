# Frozen algorithms, parameters, marker panels, and software

This document is the auditable implementation record for the integration,
label-transfer, serous-state, regulatory, and coexpression analyses. All random
operations use seed `20260831` unless explicitly stated otherwise.

## 1. Cell scores and marker panels

For cell *i* and marker panel *G*, raw counts were normalized to 10,000 counts
per cell and log transformed. The panel score was

`score(i,G) = mean[g in G]{log(1 + 10000 * count(i,g) / library_size(i))}`.

A panel required at least one detected marker. Broad epithelial assignment
required the epithelial panel to be the highest of the four broad-compartment
scores, epithelial score `>=0.15`, and at least one detected epithelial marker.
Secretory-state assignment required a winning subtype score `>=0.12` and at
least one detected winning-panel marker. If the winning-minus-runner-up margin
was `<0.05`, the cell was labelled ambiguous.

- Epithelial: `EPCAM, KRT8, KRT18, KRT19, KRT4, KRT5`
- Immune: `PTPRC, LST1, TYROBP, CD3D, CD79A, MS4A1, NKG7`
- Fibroblast: `COL1A1, COL1A2, COL3A1, DCN, LUM, COL6A1`
- Endothelial: `PECAM1, VWF, EMCN, KDR, CLDN5, RAMP2`
- Surface/club: `SCGB1A1, SCGB3A1, KRT4, KRT13, KLF3, CYP2F1`
- Serous glandular: `LTF, LYZ, BPIFB1, SLPI, WFDC2, PRR4, LPO`
- Mucous glandular: `MUC5B, BPIFB2, AGR2, SPDEF, ZG16B, PIP`
- Goblet: `MUC5AC, CLCA1, SPDEF, AGR2, FCGBP`
- Inflammatory secretory: `CST1, SLC26A4, SERPINB3, SERPINB4, CCL26, POSTN`
- Transitional secretory: `KRT4, KRT13, KRT17, KRT19, SFN, KRT16`
- Ciliary reference score: `FOXJ1, PIFO, TPPP3, CAPS, CETN2, DNAH5, DNAI2, RSPH1, CFAP54, IFT88`

All broad-compartment, epithelial, and secretory-subtype markers were placed in
the definition-gene exclusion set and were not reused as confirmatory genes in
the intrinsic-state analysis.

## 2. Integration, t-SNE, and bidirectional label transfer

The integration feature construction did not use disease labels. The common
gene universe was filtered to remove all definition genes, mitochondrial genes,
ribosomal genes, and immunoglobulin-prefix genes. A common pseudobulk variance
screen retained the top 5,000 eligible genes; required broad/subtype markers,
the ciliary panel, and the locked module were carried solely to preserve audit
fields when present in both datasets. For the actual integrated embedding,
cells were log1p-normalized to 10,000 counts per cell, and the 2,500 most
variable eligible non-definition genes were selected.

- Linear projection: randomized `TruncatedSVD`, 40 components,
  `random_state=20260831`.
- Batch correction: Harmony on the 40-dimensional projection, batch key
  `dataset`, `random_state=20260831`; other `harmonypy` arguments used package
  defaults.
- t-SNE display: corrected components 1-30; two dimensions; Barnes-Hut method;
  perplexity 50; learning rate `auto`; initialization `pca`; 1,200 iterations;
  `random_state=20260831`; `n_jobs=-1`.
- Bidirectional transfer: 25-nearest-neighbor classifier in corrected components
  1-30; distance weighting; scikit-learn default Minkowski distance (`p=2`,
  Euclidean); each dataset was used once as reference and once as query.
- Reporting: confusion matrices were row-normalized within the true/reference
  identity. Balanced accuracy, macro-F1, and serous recall were computed from
  the original versus transferred labels.

The t-SNE was used only for display. Donor-level inference did not use t-SNE
coordinates.

## 3. Serous-state coordinate

Only cells labelled serous glandular were used. Definition, mitochondrial,
ribosomal, and immunoglobulin-prefix genes were excluded. The 1,500 most
variable eligible genes were projected with randomized `TruncatedSVD` to 30
components (`random_state=20260831`) and corrected by Harmony for `dataset`.
A 25-nearest-neighbor Euclidean graph was built in corrected components 1-25.
Let `d(i,j)` denote a neighbor distance and let `sigma` be the median 25th-
neighbor distance over cells. Directed neighbor weights were

`w(i,j) = exp[-d(i,j)^2 / max(sigma^2, 1e-8)]`.

The graph was symmetrized by `(W + W^T)/2`. With degree matrix `D`, the
normalized kernel was `K = D^(-1/2) W D^(-1/2)`. The 11 largest-algebraic
eigenpairs were obtained with ARPACK `eigsh` using a start vector drawn from
`numpy.random.default_rng(20260831).normal(size=n_cells)`, `tol=1e-8`, and
`maxiter=max(1000, 5*n_cells)`. The trivial first eigenpair was discarded and
diffusion coordinate `k` was `DC_k(i) = lambda_k * psi_k(i)` for the next ten
eigenpairs.

The mature score was the arithmetic mean of log1p-normalized expression across
the detected locked-module genes. The ciliary score was the corresponding mean
for the 10-gene ciliary panel above. The root was the healthy GSE235711 serous
cell maximizing `mature_score - ciliary_score`. For cell *i*,

`distance(i) = sqrt(sum[k=1..10]{(DC_k(i) - DC_k(root))^2})`.

The reported serous-state coordinate was the average-tie percentile rank of
`distance(i)` among all serous cells. Donor/tissue inference used the median
cell coordinate. It is a descriptive diffusion-distance axis, not pseudotime or
lineage direction.

After fixing `v0`, recomputation against the archived coordinate gave cell-level
Spearman rho `0.999999999996`, maximum absolute coordinate difference
`6.78e-05`, and zero median difference across the 49 donor/tissue medians. The
two formal contrasts and their P values were unchanged.

## 4. CollecTRI and PROGENy activity scores

Within each dataset, each measured pseudobulk target gene was standardized
across eligible serous donor/tissue samples. For source or pathway *s*,

`activity(s,j) = sum[g]{z(g,j) * weight(s,g)} / sqrt(sum[g]{weight(s,g)^2})`.

- CollecTRI: duplicated source-target pairs were removed; signed frozen weights
  were retained; at least 20 measured targets were required.
- PROGENy: annotations required prior `p_value < 0.05`; within each pathway, the
  500 smallest-prior-P target genes were retained; duplicated source-target
  pairs were removed; at least 25 measured targets were required.
- Paired cohort: one-sample t test on four within-donor polyp-minus-ethmoid
  activity differences; leave-one-pair-out stability was the fraction of four
  omission estimates retaining the full-data sign.
- External cohort: Welch t test for polyp versus CRSwNP inferior turbinate.
- Priority rule: concordant direction, leave-one-pair-out sign stability
  `>=0.80`, paired `P<0.10`, external `P<0.05`, and Fisher-combined BH FDR
  `<0.10`. CollecTRI sources additionally required source-gene detection in
  `>=25%` of eligible pseudobulks in both cohorts, with detection defined as
  CPM `>=1` in the `log2(CPM+0.5)` matrix.
- FDR families: CollecTRI and PROGENy candidates were adjusted separately
  within their respective tested-activity families for the reported
  Fisher-combined FDR; the
  functional-enrichment analysis used separate library-by-direction families.

## 5. Frozen 198-gene coexpression module

Expression was residualized separately for every gene by subtracting its mean
within each anatomical-disease group. In discovery GSE136825, the residualized
locked-module score was the sample-wise mean of row-z-scored values for the 38
locked genes measured on both platforms. Every other common gene was correlated
with that score. The 160 non-module genes with the largest positive discovery
correlation were frozen, giving 198 genes total. No validation statistic was
used for neighbor selection.

- Gene membership: Pearson correlation of each residualized gene with the
  residualized standardized locked-module score, recalculated separately in
  each cohort.
- Complete pairwise topology: all `198 choose 2 = 19,503` Pearson correlations.
- Consensus edge: discovery correlation `>=0.55` and validation correlation
  `>=0.30`.
- Consensus edge weight: geometric mean
  `sqrt(r_discovery * r_validation)`.
- Density: mean absolute Pearson correlation over all 19,503 gene pairs.

The original null used 500 same-sized sets drawn without replacement from all
eligible validation genes. The strengthened sensitivity null used 5,000 draws
with seed `20260831`. Eligible genes required finite mean expression, finite
positive within-group residual variance, and exclusion from the frozen module;
the resulting non-module universe contained 18,052 genes. Genes were assigned
to a 5 x 5 joint rank-quantile grid using validation-cohort mean expression and
residual variance. Every random set reproduced the frozen module's exact gene
count in each occupied stratum and sampled without replacement within strata.
The statistic was the same validation mean absolute pairwise correlation, and
the empirical one-sided P value was

`P = (1 + number of null densities >= observed density) / (1 + number of draws)`.

The observed validation density was `0.37884065`; none of 5,000 matched random
sets reached it, giving `P=1/5001=0.00019996`.

## 6. Software and environment status

The advanced analysis environment used Python 3.12.13 with NumPy 2.4.6,
SciPy 1.15.3, pandas 3.0.5, scikit-learn 1.7.1, anndata 0.12.2,
harmonypy 0.0.10, h5py 3.16.0, statsmodels 0.14.5, networkx 3.5,
matplotlib 3.11.1, seaborn 0.13.2, and decoupler 2.1.1. Publication figures
used R 4.6.1, ggplot2 4.0.3, ragg 1.5.2, scales 1.4.0,
ComplexHeatmap 2.28.0, circlize 0.4.18, igraph 2.3.3, and ggrepel 0.9.8.

The early classification directory retained NumPy 2.5.2 and SciPy 1.18.1 but
did not retain a complete `pip freeze` or pandas distribution metadata. It is
therefore documented as an observed historical environment, not misrepresented
as a complete lock. The reproducibility folder contains:

- `requirements-advanced-full.txt`: complete advanced environment lock;
- `requirements-core-observed.txt`: explicitly incomplete historical audit;
- `environment-unified-python.yml`: proposed clean end-to-end reproduction
  environment;
- `R-packages-lock.txt`: direct R package versions used for figure production.
