# 15-Minute City Framework and Physical Activity

This repository contains the code and analysis scripts for the research paper:

> **A 15-minute city framework to active mobility and physical activity: A cross-sectional study in the Canton of Geneva, Switzerland**

## Project Overview

This study investigates the relationships between urban accessibility (operationalized through the 15-minute city framework) and physical activity patterns in Geneva, Switzerland, using data from the Bus Santé population-based study. The analysis combines comprehensive accessibility metrics with advanced spatial modeling techniques to examine both active mobility and leisure-time physical activity.

Key research questions:
1. How can we characterize and quantify spatial accessibility to diverse amenities using a comprehensive 15-minute city framework?
2. What are the associations between accessibility metrics and different domains of physical activity, particularly non-occupational active mobility and leisure-time moderate-to-vigorous physical activity?
3. Do these relationships show nonlinear patterns and threshold effects?

## Repository Structure

```
├── src/                    # Source code
│   ├── .ipynb_checkpoints/ # Jupyter notebook checkpoints
│   ├── archives/           # Archived analysis files
│   ├── cache/              # Cached computation results
│   ├── Rhistory/           # R command history
│   ├── src.Rproj           # R project file
│   ├── load_packages.R     # Package loading and environment setup
│   ├── 4. INLA Models.R    # Bayesian spatial models with SPDE approach
│   ├── make_tables_15min.R # Table generation and descriptive statistics
│   ├── utils.R             # Utility functions for analysis and visualization
│   ├── 1. OSM POI extraction and cleaning.ipynb     # OpenStreetMap POI extraction and classification
│   ├── 2. Accessibility metrics computation.ipynb   # Proximity time calculations and spatial analysis
│   └── 3. Association with physical activity.ipynb # Physical activity analysis and modeling
├── cache/                  # Cached data and intermediate results
├── data/                   # Data files (gitignored)
├── manuscript/             # Manuscript files and drafts
├── results/                # Output files and figures (gitignored)
├── .gitignore              # Git ignore specifications
├── .Rhistory               # R session history
└── README.md               # This file
```

## Analysis Workflow

The analysis follows these main steps:

1. **Data Preparation and POI Classification**:
   - OpenStreetMap data extraction and POI classification (`1. OSM POI extraction and cleaning.ipynb`)
   - POI categorization into 9 categories following established frameworks
   - Pedestrian network construction using OSMnx

2. **Accessibility Analysis**:
   - Network-based walking distance calculations (`2. Accessibility metrics computation.ipynb`)
   - Proximity time (PT) computation using Pandana
   - H3 hexagonal aggregation for spatial representation

3. **Physical Activity Analysis**:
   - PAFQ questionnaire data preparation and analysis (`3. Association with physical activity.ipynb`)
   - Non-occupational active mobility time calculation
   - Leisure-time MVPA time standardization
   - Initial statistical associations

4. **Advanced Statistical Modeling**:
   - Bayesian spatial models with SPDE approach (`4. INLA Models.R`)
   - Two-part hurdle models for zero-inflated outcomes
   - Progressive model adjustment for confounders
   - RW2 smoothing for nonlinear relationships and threshold identification

5. **Table Generation and Reporting**:
   - Descriptive statistics and publication tables (`make_tables_15min.R`)
   - Model summaries and effect estimates
   - Sensitivity analyses and robustness checks

6. **Utilities and Support**:
   - Package loading and environment setup (`load_packages.R`)
   - Utility functions for analysis and visualization (`utils.R`)
   - Cached results and intermediate computations (`cache/`)

## Key Findings

The analysis revealed:

- **Accessibility Distribution**: 59.6% of participants had 15-minute access to amenities, with clear urban-rural gradients
- **Domain-Specific Effects**: Greater accessibility increased non-occupational active mobility (2.7% per SD improvement) but decreased leisure-time MVPA (3.0% decrease per SD improvement)
- **Critical Threshold**: 30-minute proximity time marks where behavioral effects fundamentally shift from positive to negative (for active mobility) or negative to positive (for MVPA)
- **Spatial Patterns**: Three distinct accessibility patterns identified: homogeneous (transport, outdoor spaces), intermediate (supplies, services), and centralized (healthcare, education, culture)
- **Policy Implications**: 15-minute city interventions simultaneously promote active transportation while potentially constraining recreational activity

## Dependencies

The analysis uses the following main libraries:

**Python (Jupyter Notebooks):**
- osmnx, pandana (for network analysis and accessibility calculations)
- geopandas, shapely (for spatial data handling)
- h3-pandas (for hexagonal spatial indexing)
- pandas, numpy (for data manipulation)
- matplotlib, seaborn (for visualization)
- requests, json (for API interactions)

**R (Statistical Analysis):**
- INLA (for Bayesian spatial modeling with SPDE)
- sp, sf (for spatial data manipulation)
- ggplot2, viridis (for visualization)
- flextable, gtsummary (for table creation)
- dplyr, tidyr (for data manipulation)
- MASS, ordinal (for statistical modeling)

## Data Sources

The analysis integrates multiple data sources:

- **Physical Activity**: Bus Santé study (2005-2024, n=13,146 participants)
- **Points of Interest**: OpenStreetMap (9 categories, >10,000 POIs)
- **Population Data**: SITG/INSEE/SFSO population grids
- **Networks**: OSMnx pedestrian network extraction
- **Administrative**: Canton of Geneva boundaries

*Note: Individual-level data is not included due to privacy concerns.*

## Spatial Framework

The study employs a comprehensive spatial framework:

- **Study Area**: Canton of Geneva with 4km buffer to prevent edge effects
- **Accessibility Metric**: Proximity time to 20 nearest POIs per category
- **Walking Speed**: 4.5 km/h standard assumption
- **Spatial Resolution**: H3 hexagons (resolution 10, ~65m edge length)
- **Coordinate System**: Swiss LV95 (EPSG:2056)

## Contact

For questions about this code, please contact the corresponding authors.

## Citation

If you use this code in your research, please cite our paper:

```
De Ridder, D., Mechoullam, S., Lamour, J., Joost, S., Guessous, I., Specchio-COVID19 study group, & Nehme, M. (2025). A 15-minute city framework to active mobility and physical activity: A cross-sectional study in the Canton of Geneva, Switzerland.
```















