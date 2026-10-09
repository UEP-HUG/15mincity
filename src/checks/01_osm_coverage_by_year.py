"""
Check 01: OSM POI coverage over time in the Geneva study area.

Counts OSM features matching each of the nine POI category definitions used in
notebook 1, at yearly timestamps, using the ohsome API (OSM full history).
The goal is to find the year from which POI coverage plateaus, which sets the
anchor year for the time-varying exposure.

Outputs (results/checks/01_osm_coverage/):
  - counts_by_commune_year.csv : per category x commune x timestamp
  - counts_canton_year.csv     : canton totals (sum over communes)
  - counts_bbox_year.csv       : buffered bbox totals (the extent used in NB2)
  - pipeline_comparison.csv    : ohsome 2025-03-13 bbox counts vs processed parquet
  - coverage_by_category.png   : small multiples, count relative to 2025-03-13

Run from the repository root:
  ~/miniforge3/envs/15mincity/bin/python src/checks/01_osm_coverage_by_year.py

Nothing in the main pipeline is read for writing or modified.
"""
import json
import time
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import requests
from shapely.geometry import box

ROOT = Path(__file__).resolve().parents[2]
NB1 = ROOT / "src" / "1. OSM POI extraction and cleaning.ipynb"
OUT = ROOT / "results" / "checks" / "01_osm_coverage"
OUT.mkdir(parents=True, exist_ok=True)

OHSOME = "https://api.ohsome.org/v1"
CATEGORIES = ["physical", "culture", "education", "healthcare", "services",
              "transport", "outdoor", "supplies", "restaurant"]
# Mid-year snapshots (closest to the survey-year midpoint), plus the pipeline date
TIMESTAMPS = [f"{y}-07-01" for y in range(2008, 2025)] + ["2025-03-13"]
# Post-hoc exclusions applied in NB1 cell 22
EXCLUSIONS = {"physical": '"leisure"="swimming_pool"', "outdoor": '"leisure"="garden"'}


def load_tag_dicts():
    """Execute the NB1 cell that defines the *_pois_dict objects (single source of truth)."""
    nb = json.loads(NB1.read_text())
    ns = {}
    for cell in nb["cells"]:
        src = "".join(cell["source"])
        if cell["cell_type"] == "code" and "healthcare_pois_dict = {" in src:
            exec(src, ns)
            break
    return {c: ns[f"{c}_pois_dict"] for c in CATEGORIES}


def q(s):
    return '"' + str(s).replace('"', '\\"') + '"'


def ohsome_filter(tag_dict, exclusion=None):
    parts = []
    for key, values in tag_dict.items():
        if values is True:
            parts.append(f"{q(key)}=*")
        else:
            parts.append(f"{q(key)} in ({', '.join(q(v) for v in values)})")
    flt = "(" + " or ".join(parts) + ")"
    if exclusion:
        flt += f" and not {exclusion}"
    return flt


def ohsome_post(endpoint, data, retries=4):
    for attempt in range(retries):
        r = requests.post(f"{OHSOME}{endpoint}", data=data, timeout=900)
        if r.status_code == 200:
            return r.json()
        print(f"  HTTP {r.status_code}: {r.text[:300]}")
        time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"ohsome request failed: {endpoint}")


def study_geometries():
    canton = gpd.read_file(ROOT / "data/input/canton_ge.GeoJSON").to_crs(2056)
    communes = gpd.read_file(ROOT / "data/input/communes_ge.GeoJSON").to_crs(2056)
    # Same construction as NB2: bbox of the canton buffered by 4 km
    bbox = box(*canton.buffer(4000).to_crs(4326).total_bounds)
    communes = communes[["NAME", "geometry"]].dissolve("NAME").reset_index()
    communes["geometry"] = communes.simplify(20)
    communes = communes.to_crs(4326)
    fc = {"type": "FeatureCollection", "features": [
        {"type": "Feature", "id": row.NAME, "properties": {},
         "geometry": row.geometry.__geo_interface__}
        for row in communes.itertuples()]}
    return bbox, json.dumps(fc)


def main():
    tag_dicts = load_tag_dicts()
    bbox, communes_geojson = study_geometries()
    bbox_str = ",".join(f"{v:.6f}" for v in bbox.bounds)
    time_str = ",".join(TIMESTAMPS)

    rows_commune, rows_bbox = [], []
    for cat in CATEGORIES:
        flt = ohsome_filter(tag_dicts[cat], EXCLUSIONS.get(cat))
        print(f"{cat}: querying ohsome ({len(flt)} chars filter)")

        res = ohsome_post("/elements/count", {"bboxes": bbox_str, "time": time_str, "filter": flt})
        for r in res["result"]:
            rows_bbox.append({"category": cat, "timestamp": r["timestamp"][:10], "count": r["value"]})

        res = ohsome_post("/elements/count/groupBy/boundary",
                          {"bpolys": communes_geojson, "time": time_str, "filter": flt})
        for grp in res["groupByResult"]:
            for r in grp["result"]:
                rows_commune.append({"category": cat, "commune": grp["groupByObject"],
                                     "timestamp": r["timestamp"][:10], "count": r["value"]})
        time.sleep(2)

    df_c = pd.DataFrame(rows_commune)
    df_b = pd.DataFrame(rows_bbox)
    df_c.to_csv(OUT / "counts_by_commune_year.csv", index=False)
    df_b.to_csv(OUT / "counts_bbox_year.csv", index=False)
    df_k = df_c.groupby(["category", "timestamp"], as_index=False)["count"].sum()
    df_k.to_csv(OUT / "counts_canton_year.csv", index=False)

    # Sanity check: does the ohsome filter reproduce the pipeline's POI set?
    cmp = []
    bbox_2056 = gpd.GeoSeries([bbox], crs=4326).to_crs(2056).iloc[0]
    for cat in CATEGORIES:
        p = gpd.read_parquet(ROOT / f"data/output/processed_osm_{cat}_pois.parquet")
        n_pipe = int(p.to_crs(2056).within(bbox_2056).sum())
        n_ohs = int(df_b.query("category == @cat and timestamp == '2025-03-13'")["count"].iloc[0])
        cmp.append({"category": cat, "pipeline_2025_03_13": n_pipe, "ohsome_2025_03_13": n_ohs,
                    "ratio_ohsome_to_pipeline": round(n_ohs / n_pipe, 3) if n_pipe else None})
    cmp = pd.DataFrame(cmp)
    cmp.to_csv(OUT / "pipeline_comparison.csv", index=False)
    print(cmp.to_string(index=False))

    plot_coverage(df_b, df_k)


def plot_coverage(df_b, df_k):
    """Small multiples: one panel per category, count relative to the 2025-03-13 count."""
    fig, axes = plt.subplots(3, 3, figsize=(11, 9), sharex=True, sharey=True)
    for ax, cat in zip(axes.flat, CATEGORIES):
        for df, color, label in [(df_k, "#2a78d6", "Canton"), (df_b, "#eb6834", "Buffered bbox")]:
            d = df[df.category == cat].sort_values("timestamp")
            ref = d.loc[d.timestamp == "2025-03-13", "count"].iloc[0]
            years = pd.to_datetime(d.timestamp).dt.year + (pd.to_datetime(d.timestamp).dt.dayofyear / 366)
            ax.plot(years, 100 * d["count"] / ref, color=color, lw=2, label=label)
        ax.axhline(90, color="#9a9a9a", lw=0.8, ls="--")
        ax.axvline(2012.5, color="#9a9a9a", lw=0.8, ls=":")
        ax.set_title(cat.capitalize(), fontsize=10, loc="left")
        ax.set_ylim(0, 110)
        ax.grid(alpha=0.2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("% of 2025-03-13 count")
    axes[0, 0].legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("OSM POI coverage over time (dashed: 90% of 2025; dotted: mid-2012)", fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(OUT / "coverage_by_category.png", dpi=200)


if __name__ == "__main__":
    main()
