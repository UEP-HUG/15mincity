"""
Check 04: who is removed by the nearest-neighbour (< 1125 m) filter in NB3?

Re-executes the data-preparation code cells of
'src/3. Association with physical activity.ipynb' up to the spatial join
(cell 81), skipping display-only cells, then compares participants kept vs
removed by the filter applied in NB3 cells 84-85.

Outputs aggregate statistics only (no individual records) to
results/checks/04_nn_filter/.

Run from src/ (the notebook uses ../data paths):
  cd src && ~/miniforge3/envs/15mincity/bin/python checks/04_nn_filter_exclusions.py
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

SRC = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SRC))  # NB3 imports the local pyspace module
NB3 = SRC / "3. Association with physical activity.ipynb"
OUT = SRC.parent / "results" / "checks" / "04_nn_filter"
OUT.mkdir(parents=True, exist_ok=True)

LAST_CELL = 81          # pyspace.add_random_noise (jitter), just before the filter
SKIP = {79, 80}         # display-only cells (single-record print, folium map)


def run_notebook_prefix():
    nb = json.loads(NB3.read_text())
    ns = {}
    for i, cell in enumerate(nb["cells"][: LAST_CELL + 1]):
        if cell["cell_type"] != "code" or i in SKIP:
            continue
        exec("".join(cell["source"]), ns)
    return ns["gdf_pafq_final_ge_15min"]


def main():
    gdf = run_notebook_prefix()
    coords = np.column_stack([gdf.geometry.x, gdf.geometry.y])
    dist, _ = cKDTree(coords).query(coords, k=2)
    gdf["nn_dist"] = dist[:, 1]
    gdf["kept"] = gdf["nn_dist"] < 1125

    pt = "overall_15min_city_proximity_time"
    gdf["year"] = gdf["appointment_year"] + 1997
    summary = gdf.groupby("kept").agg(
        n=(pt, "size"),
        pt_mean=(pt, "mean"),
        pt_median=(pt, "median"),
        pt_p90=(pt, lambda x: x.quantile(0.9)),
        pt_max=(pt, "max"),
        nn_dist_median=("nn_dist", "median"),
    ).round(2)
    print(f"Pre-filter n = {len(gdf)}; kept = {gdf.kept.sum()} (NB3 final should be 13,146)")
    print(summary.to_string())

    bins = [0, 10, 15, 30, 1000]
    labels = ["<10", "10-15", "15-30", ">30"]
    gdf["pt_bin"] = pd.cut(gdf[pt], bins=bins, labels=labels)
    by_bin = pd.crosstab(gdf["pt_bin"], gdf["kept"], margins=True)
    by_bin["pct_removed"] = (100 * by_bin.get(False, 0) / by_bin["All"]).round(2)
    print("\nRemoved by PT bin:\n", by_bin.to_string())

    summary.to_csv(OUT / "kept_vs_removed_summary.csv")
    by_bin.to_csv(OUT / "removed_by_pt_bin.csv")


if __name__ == "__main__":
    main()
