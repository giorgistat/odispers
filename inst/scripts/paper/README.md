# Scripts for the paper analysis

These scripts reproduce the analysis and figures of *A river-network exposure model for onchocerciasis prevalence*. That analysis uses REMO nodule data for Sudan and South Sudan, with the dispersal range fixed at 40 km, breeding intensity quadratic in log discharge, and aridity in the intensity of transmission as a linear spline in log AI with a change point.

All inputs come with the package, so the scripts need no external files:

```r
data(remo)            # REMO nodule surveys, 16 countries (Sudan: country == "SUDAN")
data(sudan_rivers)    # HydroRIVERS reaches with discharge >= 0.1 m3/s (sf)
data(sudan_aridity)   # Global Aridity Index v3.1, packed terra raster (AI x 10,000)
data(sudan_borders)   # borders of Sudan and South Sudan (sf)
```

Run the scripts in order from any working folder. Each script loads the shared settings from the installed package with

```r
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
```

and writes its results to `paper_results/` and `paper_figures/` in the current working folder. To run a script without cloning the repository, use, for example, `source(system.file("scripts/paper/02_prepare_inputs.R", package = "odispers"))`.

| Script | What it does |
|---|---|
| `00_settings.R` | Paths and settings shared by the other scripts (sourced by them) |
| `01_exploratory_maps.R` | Exploratory maps and the village table (river flow, aridity, distances to rivers); standalone |
| `02_prepare_inputs.R` | Villages with aridity, quadrature nodes, village–node pairs, prediction grid |
| `03_fit_spline.R` | Main fits at 40 km and 53 km, with a profile likelihood for the change point |
| `04_river_size_and_aridity.R` | Effect of river size on breeding, share of exposure by river size, aridity curve, change-point profile |
| `05_observed_vs_predicted.R` | Observed against predicted prevalence |
| `06_prediction_maps.R` | Maps of predicted prevalence, exposure and intensity of transmission |
| `07_fit_u_field.R` | Adds a Whittle–Matérn breeding field u on the river network and a village effect to the spline/quadratic model; profiles the range of u |
| `08_fit_u_linear.R` | Simplified model with u: aridity linear in log AI and breeding linear in log discharge |
| `09_variogram_z.R` | Variogram of the predicted village effects with a permutation envelope of spatial independence, with and without u |
| `10_fit_s_vs_u.R` | Spatial random effect in the exposure (u on the river network) or in the intensity of transmission (Gaussian process S on the plane), each with and without the village effect |

Intermediate results are written to `paper_results/` and figures to `paper_figures/`. Each script uses results of earlier scripts from `paper_results/` in the working folder:

| Script | Needs |
|---|---|
| 03, 08 | `inputs.rds` (from 02; built automatically if missing) |
| 04, 05, 06, 07 | `inputs.rds` and `spline_fits.rds` (from 03) |
| 09 | `u_fits_linear.rds` (from 08) |
| 10 | `inputs.rds`; uses `u_fits_linear.rds` from 08 if present |

So the linear model with u can be run on its own with `08_fit_u_linear.R`. Script 03 fits the model about 110 times to profile the change point and takes some minutes; scripts 07 and 08 take about an hour each.

`../extract_riveratlas_sudan.R` extracts RiverATLAS attributes for the study area, such as stream gradient and minimum and maximum monthly discharge. These are for planned breeding covariates.
