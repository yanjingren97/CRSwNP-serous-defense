# Reproducibility environments

The project contains two historically distinct Python dependency directories.
They are documented separately rather than being retrospectively represented as
one environment.

- `requirements-core-observed.txt` records the distribution metadata still
  present for the early classification environment. It is incomplete and is an
  audit record, not a complete lock file.
- `requirements-advanced-full.txt` records every distribution in the advanced
  analysis environment used for integration, diffusion-state, regulatory, and
  coexpression analyses.
- `environment-unified-python.yml` is the proposed clean Python environment for
  future end-to-end reproduction. It uses the advanced analysis versions and
  must be identified as a reproduction environment rather than the historical
  environment.
- `R-packages-lock.txt` records the packages directly used to render the main
  figures. Full transitive R session information remains in each figure package.

Python version used for the advanced audit: 3.12.13, MSC v.1944, 64 bit.
R version used for publication figures: 4.6.1 on Windows 11 x64.
