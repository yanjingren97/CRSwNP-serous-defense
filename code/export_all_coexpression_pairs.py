from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from bulk_coexpression_preservation import load_gse136825, load_gse36830, residualize


def main():
    out_dir = ROOT / "results" / "advanced_coexpression"
    nodes = pd.read_csv(out_dir / "preserved_module_nodes.csv")
    module = nodes["gene"].tolist()

    discovery_expr, discovery_group = load_gse136825()
    validation_expr, validation_group = load_gse36830()
    common = discovery_expr.index.intersection(validation_expr.index)
    discovery_expr = discovery_expr.loc[common]
    validation_expr = validation_expr.loc[common]

    missing = sorted(set(module).difference(common))
    if missing:
        raise RuntimeError(f"Module genes missing from common expression matrix: {missing}")

    discovery_residual = residualize(discovery_expr, discovery_group)
    validation_residual = residualize(validation_expr, validation_group)
    discovery_cor = np.corrcoef(discovery_residual.loc[module].to_numpy())
    validation_cor = np.corrcoef(validation_residual.loc[module].to_numpy())
    tri = np.triu_indices(len(module), 1)

    pairs = pd.DataFrame({
        "source": np.asarray(module, dtype=object)[tri[0]],
        "target": np.asarray(module, dtype=object)[tri[1]],
        "cor_discovery": discovery_cor[tri],
        "cor_validation": validation_cor[tri],
    })
    pairs["is_consensus"] = (
        (pairs.cor_discovery >= 0.55) &
        (pairs.cor_validation >= 0.30)
    )
    pairs["consensus_weight"] = np.nan
    consensus_mask = pairs.is_consensus
    pairs.loc[consensus_mask, "consensus_weight"] = np.sqrt(
        pairs.loc[consensus_mask, "cor_discovery"] *
        pairs.loc[consensus_mask, "cor_validation"]
    )

    rho, p_value = spearmanr(
        pairs.cor_discovery,
        pairs.cor_validation,
        nan_policy="omit",
    )
    if len(pairs) != 19_503:
        raise RuntimeError(f"Expected 19,503 pairs, found {len(pairs):,}")
    if int(pairs.is_consensus.sum()) != 7_389:
        raise RuntimeError(
            f"Expected 7,389 consensus edges, found {int(pairs.is_consensus.sum()):,}"
        )

    output = out_dir / "all_module_gene_pairs.csv"
    pairs.to_csv(output, index=False)
    print(
        f"pairs={len(pairs):,}; consensus={int(pairs.is_consensus.sum()):,}; "
        f"Spearman rho={rho:.15f}; p={p_value:.15e}"
    )
    print(output)


if __name__ == "__main__":
    main()
