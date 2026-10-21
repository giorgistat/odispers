# Scripts for the paper analysis

These scripts reproduce the analysis and figures of *A river-network exposure model for onchocerciasis prevalence*. That analysis uses REMO nodule data for Sudan and South Sudan, with the dispersal range fixed at 40 km, breeding intensity quadratic in log discharge, and aridity in the intensity of transmission as a linear spline in log AI with a change point.

Run them in order from a folder that contains:
- the REMO data, `Remo_AllCountries.csv`;
- HydroRIVERS for Africa, `HydroRIVERS_v10_af_shp/`;
- the annual Global Aridity Index v3.1 GeoTIFF.

| Script | What it does |
|---|---|
| `00_settings.R` | Paths and settings shared by the other scripts (sourced by them) |
| `01_exploratory_maps.R` | Exploratory maps and the village table (river flow, aridity, distances to rivers); standalone |
| `02_prepare_inputs.R` | Villages with aridity, quadrature nodes, village–node pairs, prediction grid |
| `03_fit_spline.R` | Main fits at 40 km and 53 km, with a profile likelihood for the change point |
| `04_river_size_and_aridity.R` | Effect of river size on breeding, share of exposure by river size, aridity curve, change-point profile |
| `05_observed_vs_predicted.R` | Observed against predicted prevalence |
| `06_prediction_maps.R` | Maps of predicted prevalence, exposure and intensity of transmission |

Intermediate results are written to `paper_results/` and figures to `paper_figures/`. Script 03 fits the model about 110 times to profile the change point and takes some minutes.

`../extract_riveratlas_sudan.R` extracts RiverATLAS attributes for the study area, such as stream gradient and minimum and maximum monthly discharge. These are for planned breeding covariates.
