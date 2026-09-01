from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def main():
    base = ROOT / "results" / "advanced_singlecell"
    old = pd.read_csv(base / "pre_deterministic_v0_fix" / "serous_cell_state_axis.csv").set_index("cell_id")
    new = pd.read_csv(base / "serous_cell_state_axis.csv").set_index("cell_id")
    cells = old[["failure_axis"]].join(new[["failure_axis"]], lsuffix="_old", rsuffix="_new", how="inner")

    old_donor = pd.read_csv(base / "pre_deterministic_v0_fix" / "serous_donor_state_axis.csv")
    new_donor = pd.read_csv(base / "serous_donor_state_axis.csv")
    keys = ["dataset", "donor", "disease", "tissue"]
    donors = old_donor.merge(new_donor, on=keys, suffixes=("_old", "_new"))
    cell_delta = np.abs(cells.failure_axis_old - cells.failure_axis_new)
    donor_delta = np.abs(donors.median_failure_axis_old - donors.median_failure_axis_new)

    audit = pd.DataFrame(
        [
            {"metric": "n_cells_compared", "value": len(cells)},
            {"metric": "cell_failure_axis_spearman_rho", "value": spearmanr(cells.failure_axis_old, cells.failure_axis_new).statistic},
            {"metric": "cell_failure_axis_max_absolute_difference", "value": cell_delta.max()},
            {"metric": "cell_failure_axis_median_absolute_difference", "value": np.median(cell_delta)},
            {"metric": "n_donor_tissue_units_compared", "value": len(donors)},
            {"metric": "donor_median_max_absolute_difference", "value": donor_delta.max()},
            {"metric": "donor_median_median_absolute_difference", "value": np.median(donor_delta)},
            {"metric": "eigensolver_seed", "value": 20260831},
            {"metric": "eigensolver_tol", "value": 1e-8},
        ]
    )
    audit.to_csv(base / "deterministic_v0_reproducibility_audit.csv", index=False)
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
