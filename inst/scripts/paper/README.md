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
| `11_fit_u_s_z.R` | Full model with u in the exposure, S in the intensity and the village effect |
| `12_cv_calibration.R` | Cross-validated calibration of the model of 08 (10 random folds): non-randomised PIT among villages with at least one positive, and discrimination and reliability of the probability of at least one positive |
| `13_fit_zero_inflated.R` | Zero-inflated model: probability of a transmission zone 1 - exp(-alpha E) and force of infection E^gamma T, with u and the village effect; gamma free, 0 and 1 |
| `14_cv_zero_inflated.R` | Cross-validated calibration of the zero-inflated model, compared with the model without zero inflation |
| `15_fit_zi_aridity.R` | Zero-inflated model with aridity also in the probability of a transmission zone, 1 - exp(-alpha E AI^delta), and its cross-validated calibration |
| `16_zi_report_figures.R` | Report figures for the zero-inflated model: maps of the probability of a transmission zone and of expected prevalence, and random-fold calibration |
| `17_fit_zil.R` | Zero-inflated model with log river exposure and log aridity in both parts (logit for the zone of transmission, log for the force of infection), on the full river network; maps and random-fold calibration |
| `18_variogram_z_zil.R` | Variogram of the village effect in the model of 17 among villages with at least one positive, with a 9999-permutation envelope, a global test and a one-sided test of positive correlation; observed against predicted prevalence, in-sample and cross-validated |
| `19_zil_report.R` | Figures for the stand-alone report on the model of 17 (`inst/reports/zero_inflated/`): maps of u, exposure on the grid and along the river network, probability of a zone, prevalence within a zone and expected prevalence, and the cross-validated non-randomised PIT and reliability of this model alone |
| `21_fit_directed.R` | Directed breeding fields (tail-down, tail-up with discharge-weighted mixing, and their sum; `directed_graph()`) in the model of 17, compared with the symmetric field over ranges of u, with random-fold cross-validation (several hours; restartable) |
| `22_profile_phi.R` | Profile of the likelihood over the tail-up share of variance in the tail-up + tail-down field, at the tail-down estimates (needs `directed_taildown.rds` from 21) |
| `23_directed_profiles.R` | Fair comparison of the directed fields (symmetric field dropped): dense range profiles, joint estimation of the range (separate ranges for the two parts of the summed field), random-fold cross-validation and comparison; restartable from checkpoints |
| `24_tailup_joint_warm.R` | Joint fit of the tail-up model started from the latent mode of its best profile fit (its latent posterior has several modes) |
| `25_directed_maps.R` | Maps of u and of expected prevalence for the joint fits of the three directed fields |
| `20_zil_intervals.R` | 95% confidence intervals for all parameters of the model of 17: Wald intervals for the regression coefficients, and for a_E, gamma, sigma_u and sigma_Z from the finite-difference observed information of the approximate marginal likelihood (about 15 minutes) |

Intermediate results are written to `paper_results/` and figures to `paper_figures/`. Each script uses results of earlier scripts from `paper_results/` in the working folder:

| Script | Needs |
|---|---|
| 03, 08 | `inputs.rds` (from 02; built automatically if missing) |
| 04, 05, 06, 07 | `inputs.rds` and `spline_fits.rds` (from 03) |
| 09 | `u_fits_linear.rds` (from 08) |
| 10 | `inputs.rds`; uses `u_fits_linear.rds` from 08 if present |
| 11 | `inputs.rds`; compares with `s_vs_u_fits.rds` from 10 if present |
| 12 | `u_fits_linear.rds` (from 08) |
| 13 | `inputs.rds`; compares with `u_fits_linear.rds` from 08 if present |
| 14 | `zi_fits.rds` (from 13) and `cv_calibration.rds` (from 12) |
| 15 | `inputs.rds`, `cv_calibration.rds` (from 12); starts from `zi_fits.rds` (13) if present |
| 16 | `zi_fits.rds` (13) and `cv_zero_inflated.rds` (14) |
| 17 | `inputs.rds`, `cv_calibration.rds` (12); compares with `cv_zero_inflated.rds` (14) if present |
| 18 | `zil_fit.rds` and `cv_zil.rds` (17) |
| 19 | `zil_fit.rds`, `cv_zil.rds` (17) and `inputs_full_network.rds` (17) |
| 20 | `zil_fit.rds` and `inputs_full_network.rds` (17) |
| 21 | `zil_fit.rds`, `cv_zil.rds`, `inputs_full_network.rds` (17) and `cv_calibration.rds` (12) |

So the linear model with u can be run on its own with `08_fit_u_linear.R`. Script 03 fits the model about 110 times to profile the change point and takes some minutes; scripts 07 and 08 take about an hour each.

`../extract_riveratlas_sudan.R` extracts RiverATLAS attributes for the study area, such as stream gradient and minimum and maximum monthly discharge. These are for planned breeding covariates.
