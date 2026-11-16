# Inference for hidden river breeding sites

`odisp_fit_lgcp()` fits one joint point pattern to the village observations. Only survey outcomes and positions are observed. A site may contribute to several villages. The default model gives each site unit strength and its position is approximated by its river-cell midpoint.

The village model accepts `link = "logit"` or the default `link = "cloglog"`. Both use $\log(E_i)+x_i^\top\beta$ as the predictor, where $E_i$ sums the marked dispersal kernels. Under logit, doubling exposure doubles infection odds, and $E_i=0$ gives zero risk. The Cox source intensity keeps its log link. Add `--link=logit` to `full_river_chain.R` and use a separate output directory when comparing links. Simulations read the link from the saved fit. `Rscript check_two_point.R logit` compares sampled site counts and mark spread to quadrature.

## Village occurrence

An optional `occurrence` formula assigns one latent transmission state to each village. With \(q_i=\mathrm{logit}^{-1}(w_i^\top\gamma)\), an inactive village returns zero positives. An active village has binomial risk \(p_i\) from the marked river sites and the selected link. Its marginal probability of zero positives is \(1-q_i+q_i(1-p_i)^{n_i}\). A positive count fixes the village state to active. For a zero count, the sampler uses the binomial probability of seeing zero in an active village. The archived draws include the state when computing prevalence. To compare models, run `zero_prepare.R` with the full coordinate inputs and an output directory. `zero_chain.R` fits a model/fold/seed combination and `zero_predict.R` generates new full-river landscapes. `zero_summary.R` pools the five spatial holdouts. `zero_check.R` verifies the occurrence and site posteriors in an example with direct integration.

## Field updates

With `occurrence = ~ climate, occurrence_exposure = TRUE`, the activity logit includes a fitted coefficient times the realized raw kernel exposure. Each count, site-relocation and mark-spread proposal reevaluates the activity likelihood at its proposed exposure. Simulations use the newly generated sites in both village predictors. With zero exposure the active-village infection risk is zero for every value of the activity probability.

The integer counts and the Gaussian graph field have different conditional distributions. Given the counts, the binomial observations add no further information about the field. The graph field and any breeding coefficients have a Gaussian prior and a Poisson log likelihood. The default field update uses a Gaussian proposal at the conditional mode with a Metropolis correction. Setting `field_update = "slice"` instead draws an ellipse using the Gaussian prior and samples along it using the conditional likelihood. Both updates target the same posterior on the chosen mesh. The count updates in both fits retain the joint likelihood of all villages. [Elliptical slice sampling](https://homepages.inf.ed.ac.uk/imurray2/pub/10ess/)

## Where variational Bayes fits

A positive survey rules out every site configuration in which all cells in reach of that village are empty. With zero background exposure, a product of independent finite-rate Poisson variational factors gives that impossible configuration positive mass and hence an expected log likelihood of minus infinity.

The `check_one_cell.R` benchmark uses a distribution over entire feasible count configurations. For each count it optimizes a Gaussian approximation to the conditional field, then weights the count by the exponential of the optimized evidence bound times the village likelihood. The small example can enumerate all relevant counts. Extending this to many cells would require a way to explore feasible configurations and approximate their relative weights. The current larger-graph fitter samples these configurations. [Gaussian variational approximation for Poisson data](https://arxiv.org/abs/1709.05885)

It is also possible to sum out the sites *before* integrating the field. Conditional on intensity, the Poisson probability-generating functional gives the joint survival probability for any specified set of people. Expanding the binomial factors recovers the full likelihood. The direct expansion has up to the product of (positives + 1) over villages signed terms, with cancellation. The benchmark verifies both orders of integration for two villages.

## Reproduce the checks

Install the current odispers source, then run the following from this directory. The REMO pilot additionally uses the R package `posterior`. Each pilot scenario writes separate fits, per-chain diagnostics and village summaries.

```sh
Rscript check_one_cell.R
Rscript remo_pilot.R results reference
Rscript remo_pilot.R results slice
Rscript remo_pilot.R results half_rate
Rscript remo_pilot.R results fine_cells
```

The first script compares count means, standard deviations and probabilities of one site against direct Gaussian quadrature. The REMO script selects eight geographically neighboring surveys before inspecting their outcomes. It varies the field update, baseline source density and river-cell spacing one at a time. These checks assess implementation and sensitivity on a subnetwork. The fitted cell probabilities describe uncertainty about unmeasured sites.

## Finite marks on the full river

With `marks = "two_point"`, a site independently draws strength $1-d$ or $1+d$ with probability one half each. The fitted spread $d$ has a Uniform(0,1) prior. Conditional on $N_j$ sites in cell $j$, the number of strong sites $H_j$ is binomial with probability one half. The cell contributes $K_{ij}[(1-d)(N_j-H_j)+(1+d)H_j]$ to village $i$. A single mark applies to every village reached by a site. The field likelihood depends only on $N_j$. Given the type counts, a one-dimensional slice update estimates $d$ from the village likelihood. When $d$ is close to one, about half the sites effectively contribute very little. The inferred contrast is therefore sensitive to the assumed reference site density. [Marking theorem, chapter 5](https://stoch.math.kit.edu/english/297.php)

The count updates use the sparse kernel to find affected villages and recompute exact, nonnegative exposure sums for those villages. Optional relocation proposals keep a site's type and move it within a spatial tile. The proposal selects cells in proportion to their Poisson mean, which cancels the prior in the acceptance ratio. The tiles change during the chain. On undirected graphs, `field_update = "block_slice"` conditions on alternate bands of graph vertices and draws a separate slice angle for each independent component. The sampler shifts the bands between sweeps. Both changes preserve the discretized posterior.

### Reproduce the full runs

The scripts accept the full-data input RDS from the Sudan analysis. This minimal equivalent input can also be rebuilt from the shipped data for the intercept-only analyses.

```r
library(odispers)
library(sf)
data(remo)
data(sudan_rivers)
v <- subset(remo, country == "SUDAN")
v$y <- round(v$y)
v <- st_as_sf(v, coords = c("long", "lat"), crs = 4326, remove = FALSE)
nodes <- river_nodes(sudan_rivers[sudan_rivers$DIS_AV_CMS >= 0.1, ], spacing = 2)
saveRDS(list(villages = v, nodes = nodes, graph = river_graph(nodes, vertex_spacing = 5),
             pairs = river_pairs(v, nodes, 40)), "full_inputs.rds")
```

Prepare the standardized covariates first. Run the following from the script directory. Each chain has its own seed and output file. The `--relocate` run uses 8,000 sweeps, discards 1,500 and proposes single-cell count changes plus site relocations. Omitting the flag runs the initial 5,000-sweep comparison without relocations. The source density, field and village model are the same in both passes.

```sh
Rscript full_river_covariates.R full_inputs.rds adjusted
for model in unit two_point; do
  for seed in 8106 8206 8306 8406; do
    Rscript full_river_chain.R adjusted/inputs.rds results/observed "$model" "$seed" --relocate
  done
  Rscript full_river_summary.R adjusted/inputs.rds results/observed "$model"
done
Rscript full_river_null.R adjusted/inputs.rds results/observed results/null
Rscript full_river_null.R adjusted/inputs.rds results/observed results/marked_predictive two_point
for seed in 8106 8206 8306 8406; do
  Rscript full_river_chain.R adjusted/inputs.rds results/null_refit two_point "$seed" results/null/first_y.rds --relocate
done
Rscript full_river_summary.R adjusted/inputs.rds results/null_refit two_point
Rscript full_river_figures.R adjusted/inputs.rds results
Rscript full_river_effects.R results
Rscript check_two_point.R
```

The null redraws all site counts, including cells beyond reach of the surveys. It uses the actual survey positions and trial counts. Its simulated surveys assess the marginal spatial pattern under a new realization of the model. Village fitted means from the chains instead condition on the observed outcomes and inferred local sites. The seeded first replicate was chosen in advance for the map and null refit. The local exposure illustration holds that replicate's sites fixed while varying their strengths.

### River covariates for the primary analysis

Run `full_river_covariates.R full_inputs.rds adjusted` with Rscript before fitting. Use `adjusted/inputs.rds` for **all** subsequent chain, simulation, refit, summary and plotting calls. The same chain script selects the formulas from the prepared inputs. To reproduce the earlier intercept-only comparison, pass the original unadjusted inputs and the explicit `--intercept-only` flag to the chain script. The default full-data analysis requires covariates. With `--extended`, the runner takes 14,000 sweeps, discards 2,500 and retains site relocation moves.

The LGCP log intensity uses log discharge, log-transformed river slope, and log-transformed inundation capped at 1.5, standardized over **all** river cells with cell length as the weight. The transmission model uses standardized log village aridity. The saved scaling table defines the units of the effects and makes it possible to reproduce the intensity away from observed villages. The baseline is fixed at 0.04 sites/km at zero field and mean covariates, with no extra free source intercept. Coefficients have independent N(0, 1.5²) priors in these runs.

Conditional on site counts and the river field, the breeding coefficients form a small Poisson regression. The block sampler uses a joint Gaussian proposal based on its score and curvature and includes the reverse proposal density in its Metropolis correction. For multivariate transmission regression the block sampler mixes elliptical slices with corrected local Gaussian proposals based on the complementary log-log score and observed curvature. Each forward simulation draws the source and transmission coefficients together, evaluates the source predictor on all cells, and applies aridity at the villages before the complementary log-log survey link.

To fit a geographic trend on the river intensity, pass `--coordinates` to `full_river_covariates.R` with a separate output directory. It calculates WGS84 longitude and latitude of all river-cell midpoints and length-standardizes them before adding them to the breeding formula. Continue with those inputs for the unit and two-point fits, the simulations and the marked refit of the first simulated null survey. The village transmission formula still contains log aridity. Use `full_river_effects.R` for the coefficient plot and aridity-group predictive checks.

### Dispersal radius and graph-field range

These parameters enter different layers. At planar distance $r$ km a site contributes its strength times $3(1-r/\phi)_+^4(1+4r/\phi)/(2\phi)$. The dispersal radius $\phi$ determines which river cells reach a village. The graph-field range controls dependence of log source intensity along the river network.

Use `full_river_covariates.R full_inputs.rds out --coordinates --radius=80` with Rscript to prepare a radius-specific design. The script rebuilds the village–river pairs when the requested radius exceeds the cached search distance. Select a graph-field range independently using the chain option `--field-range=100`. For the grid, `--grid` uses 6,000 full-graph sweeps, discards 1,000 and includes local site relocations. Example for one chain:

```sh
Rscript full_river_chain.R out/inputs.rds grid_case/observed two_point 8106 --grid --field-range=100
```

The simulator reads the actual kernel radius, field range, scale and baseline rate from the fits and checks that the predictors and pair search agree. The grid compares radii 40, 60 and 80 km and field ranges 10, 25 and 100 km. The fixed 40-km distance in the neighboring-survey diagnostic is common to all models. Report site counts both in each radius's footprint and in the common 40-km footprint. The wider-kernel count updates maintain balanced trees of nonnegative contributions per village, making a changed site's exposure update logarithmic in village neighborhood size.

## Code map, regression checks and report

The sampler is organized by responsibility:

- `R/lgcp.R`: public controls, joint sampling order and retained summaries.
- `R/lgcp_setup.R`: fixed survey/river design, kernel support and prior precision.
- `R/lgcp_sampler.R`: feasible initial counts, conditional field and regression updates.
- `R/lgcp_blocks.R`: conditional graph partitions and blocked field updates.
- `R/lgcp_kernel.R`: stable truncated Gaussian normalization and width proposals.
- `src/lgcp.cpp`: local count/relocation updates and nonnegative exposure trees.

The [LGCP report](../../reports/lgcp/lgcp_river_model.pdf) follows the existing
river reports and includes the model, maps, spatial holdout scores, chain
diagnostics and the 10 November code audit. Scientific tables summarize the
archived 7 November runs; all model/fold combinations have convergence flags.
The report does not claim these are converged posterior comparisons.

From the repository root, run the regression suite and rebuild the figures:

```sh
Rscript -e 'testthat::test_local(filter = "lgcp")'
Rscript inst/scripts/lgcp/report_figures.R
```

With the current package installed, the independent integration checks are:

```sh
Rscript inst/scripts/lgcp/check_one_cell.R
Rscript inst/scripts/lgcp/check_two_point.R logit
Rscript inst/scripts/lgcp/zero_check.R
Rscript inst/scripts/lgcp/check_width.R
```

`benchmark.R` accepts a full prepared input, a directory containing an
installed version of the package, and an output prefix. Compare versions in
separate R processes so each uses its matching compiled code. The fixed
benchmark seed and settings are visible at the top of the script. It saves
four timings, a retained fit for numerical replay, and the R session details.
The 900-village benchmark decreased from 12.677 to 11.431 seconds (medians,
300 sweeps including initialization). This is a workload-specific timing,
not a measurement of effective independent posterior draws per second.

## Direct field-mixing diagnostics

The ordinary fit stores coefficient trajectories and running cell summaries,
but does not retain the Gaussian field trajectory. `field_mixing.R` records
160 selected vertex values and five spatial field summaries in four additional
tail-down, 800 km chains. It uses a private wrapper around the field update;
an identical-seed replay checks that recording changes neither the returned
fit nor the final random-number state. The existing sampler interface is unchanged.

From the repository root, with prepared inputs and the current package installed:

```sh
for chain in 1 2 3 4; do
  Rscript inst/scripts/lgcp/field_mixing.R inputs.rds library_dir field_results "$chain"
done
Rscript inst/scripts/lgcp/field_mixing_summary.R field_results
```

Each chain has 8,000 sweeps, with the first 2,000 discarded. Initial widths
are 10, 20, 40 and 80 km; every field starts at zero. Numerical diagnostics
and autocorrelations use all retained states. The trace plot alone is thinned.
The report bundles `field_diagnostics.csv`, `field_group_summary.csv` and
`field_plot_data.rds` from the 10 November run, separately from the older
forecast comparisons. `report_figures.R` also rebuilds its field figure and table.
