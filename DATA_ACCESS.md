# Data access and expected local layout

This repository does not redistribute raw GEO, GeoMx, or ProteomeXchange data. Download the public source files from their primary repositories and keep them in the ignored `data/raw/` directory.

Expected local layout for full reanalysis:

```text
data/
├── raw/
│   ├── single_cell/
│   ├── bulk_crs/
│   │   ├── GSE36830_series_matrix.txt.gz
│   │   ├── GSE136825_series_matrix.txt.gz
│   │   ├── GSE136825_genecounts_20190903.txt.gz
│   │   ├── GPL570.annot.gz
│   │   └── ...
│   ├── spatial/
│   │   └── Validation_cohort_Q3NormalizationFile.xlsx
│   └── proteomics/
│       └── PXD013330_MQ_results/proteinGroups.txt
└── resources/
    ├── gene_sets/
    └── regulatory/
        ├── CollecTRI_regulons.csv
        └── PROGENy_annotations.tsv
```

The initial single-cell preparation script retains the original extraction logic and expects downloaded supplementary archives under `data/derived/extracted/`. Because GEO supplementary filenames can change when manually downloaded or mirrored, verify filenames against `code/prepare_and_qc.py` before running.

The repository already includes frozen, non-identifying processed outputs needed to audit the main results and regenerate publication figures. Full raw-data reprocessing is therefore optional for figure recreation but required for an independent end-to-end rerun.

No private clinical data are required or included.
