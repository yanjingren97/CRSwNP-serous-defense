from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]

import anndata as ad
import numpy as np
import pandas as pd

from advanced_singlecell_integration import DEFINITION, serous_continuum


def main():
    source = ROOT / "results" / "advanced_singlecell" / "integrated_secretory_balanced.h5ad"
    adata = ad.read_h5ad(source)
    # The historical h5ad preserved the normalized sparse matrix but lost its
    # obs/var annotations during a pandas/anndata string-serialization workaround.
    # Restore the frozen row and column metadata without changing X.
    obs = pd.read_csv(
        ROOT / "results" / "advanced_singlecell" / "integrated_secretory_embedding.csv"
    ).set_index("cell_id")
    panel = pd.read_csv(
        ROOT / "results" / "advanced_singlecell" / "integration_gene_panel.csv"
    )["gene"].astype(str)
    if adata.n_obs != len(obs) or adata.n_vars != len(panel):
        raise ValueError("Frozen h5ad dimensions do not match the audit tables")
    adata.obs_names = obs.index.astype(str)
    adata.obs = obs
    adata.var_names = panel.to_numpy()
    adata.var["eligible_embedding"] = np.asarray(
        [
            g not in DEFINITION and not g.startswith(("MT-", "RPL", "RPS", "IG"))
            for g in adata.var_names
        ],
        dtype=bool,
    )
    serous_continuum(adata)
    print(f"recomputed_serous_cells={(adata.obs['subtype'].astype(str) == 'Serous_glandular').sum()}")


if __name__ == "__main__":
    main()
