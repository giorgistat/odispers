# odispers

River-network dispersal models for onchocerciasis prevalence.

`odispers` fits exposure models for prevalence survey data in which vector breeding takes place on a river network, exposure at each location accumulates from all river segments within reach, and the distance over which rivers spread infection depends on aridity. Prevalence is exactly zero at locations beyond the reach of any river, so areas without transmission follow from the model rather than being specified in advance.

The basic model has no latent fields and is fitted by maximum likelihood (`odisp_fit()`). Version 0.3 adds a Gaussian Whittle–Matérn field $u$ on the river network in the breeding intensity, together with an independent village effect, fitted by a Laplace approximation (`odisp_fit_u()`; see below). A Gaussian process on the plane will follow.

## Model

For a survey location $x$ with $n$ people examined and $y$ found positive,

$$
w(s) = \exp\{z(s)^\top\theta\}, \qquad
E(x) = \int_\Gamma K_{\phi(x)}(\lVert x - s\rVert)\, w(s)\, ds,
$$

$$
\log\phi(x) = \eta_0 + \gamma \log \mathrm{AI}(x) \quad \text{(optionally bounded above by } \phi_{\max}\text{)}, \qquad
Y \sim \mathrm{Bin}\Big(n,\; 1 - \exp\big[-\{E(x) + \epsilon\}\, e^{d(x)^\top\beta}\big]\Big),
$$

where

- $\Gamma$ is the river network and $z(s)$ are characteristics of the river at $s$, such as log discharge;
- $K_\phi(r) = k(r/\phi)/(c\,\phi^a)$, where $k(t) = (1-t)_+^4(1+4t)$ is the Wendland kernel. The kernel is zero beyond $t = 1$, so a river contributes to $E(x)$ only within the dispersal range $\phi(x)$. Two normalisations are available:
  - `scaling = "line"` ($a = 1$, $c = 2/3$): exposure on the bank of an isolated straight river equals its breeding intensity, whatever $\phi$.
  - `scaling = "area"` ($a = 2$, $c = \pi/7$): exposure is a kernel-weighted density of breeding habitat, in km of river per km$^2$ weighted by intensity, and does not grow with $\phi$ where the network is dense;
- $\mathrm{AI}(x)$ is the aridity index (precipitation divided by reference evapotranspiration), so $\gamma > 0$ means a shorter range where it is drier;
- $d(x)$ are covariates on the plane, and $\epsilon \ge 0$ is an optional small background exposure.

The integral is computed by quadrature over short pieces of the network.

## Breeding field on the river network

`odisp_fit_u()` keeps the dispersal range fixed and adds two latent effects:
- a field $u$ on the network, which enters the breeding intensity as $w(s) = \exp\{z(s)^\top\theta + u(s)\}$;
- a village effect $Z$, which enters the force of infection.

The field $u$ is a Whittle–Matérn field on the metric graph with $\alpha = 1$ (Bolin, Simas & Wallin, 2024). Its covariance is exponential along the rivers, and its precision at the graph vertices is exact and sparse (`wm_precision()`). `river_graph()` builds the graph from HydroRIVERS: reaches are joined through `NEXT_DOWN`, and reaches longer than `vertex_spacing` are split into shorter edges.

```r
graph <- river_graph(nodes, vertex_spacing = 5)        # nodes need HYRIV_ID and NEXT_DOWN
fits  <- odisp_fit_u(y ~ log(AI) + h, ~ log(DIS_AV_CMS) + I(log(DIS_AV_CMS)^2),
                     villages, nodes, graph, phi = 40, range_u = c(25, 50, 100, 200))
attr(fits, "profile")                                  # Laplace log-likelihood by range of u
```

How the fit works:
- **Integrating out the latent effects.** The regression coefficients get diffuse Gaussian priors, and they are integrated out together with $u$ and $Z$.
- **Finding the mode.** The mode is found by truncated Newton with the exact Hessian. Conjugate gradients are preconditioned with the Gauss–Newton Hessian, which is inverted through the Woodbury identity using one sparse Cholesky factor per range. For the Sudan data (75,000 vertices, 900 villages), one Laplace evaluation takes 10–20 seconds.
- **Hyperparameters.** $\sigma_u$ and $\sigma_Z$ are chosen by maximising the approximate marginal likelihood, and the range of $u$ is profiled.

The full analysis is in `inst/scripts/paper/07_fit_u_field.R`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("giorgistat/odispers")
```

The likelihood is computed in compiled C++ code (via Rcpp), so installing from GitHub needs a working compiler: Rtools on Windows, the Xcode command line tools on macOS.

## Usage

```r
library(odispers)
library(sf)

# survey locations: an sf object with columns y (positives), n (examined) and AI (aridity)
villages$AI <- extract_aridity(villages, "ai_v31_yr.tif")

# river network, e.g. HydroRIVERS, keeping reaches with some flow
rivers <- rivers[rivers$DIS_AV_CMS >= 0.1, ]

nodes <- river_nodes(rivers, spacing = 2)               # quadrature nodes every 2 km
pairs <- river_pairs(villages, nodes, max_range = 100)  # village-node pairs within 100 km

fit <- odisp_fit(y ~ 1, breeding = ~ log(DIS_AV_CMS),
                 data = villages, nodes = nodes, pairs = pairs, max_range = 100)
summary(fit)

# does the range depend on aridity?
fit0 <- update(fit, fixed = list(gamma = 0))
AIC(fit0, fit)
odisp_profile(fit, "gamma", seq(0, 1, by = 0.1))

# prediction on a grid (needs AI at each grid point)
pred <- predict(fit, grid, nodes)   # phi, E, prevalence, zero_exposure
```

A full analysis of REMO nodule data for Sudan and South Sudan, with HydroRIVERS and the Global Aridity Index, is in `inst/scripts/fit_sudan.R`.

## Data

The package includes the data for the Sudan and South Sudan analysis:

| Dataset | Content | Source |
|---|---|---|
| `remo` | REMO nodule surveys, 14,126 villages in 16 countries | APOC REMO (Zouré et al., 2014) |
| `sudan_rivers` | 93,094 HydroRIVERS reaches with discharge ≥ 0.1 m³/s, with `NEXT_DOWN` topology | Lehner & Grill (2013) |
| `sudan_aridity` | Annual aridity index at 30″, as a packed `terra` raster (values × 10,000) | Global Aridity Index v3.1, Zomer et al. (2022), CC BY 4.0 |
| `sudan_borders` | Borders of Sudan and South Sudan | Natural Earth |

```r
data(remo); data(sudan_rivers); data(sudan_aridity)
villages <- sf::st_as_sf(subset(remo, country == "SUDAN"), coords = c("long", "lat"), crs = 4326)
villages$AI <- extract_aridity(villages, sudan_aridity)
```

`data-raw/make_data.R` rebuilds the datasets from the original files.

## Main functions

| Function | Purpose |
|---|---|
| `river_nodes()` | Splits the river network into quadrature nodes with their lengths and attributes |
| `river_pairs()` | Finds the nodes within `max_range` of each location |
| `aggregate_nodes()` | Merges nodes onto a coarser grid, within classes of flow, to speed up fitting |
| `extract_aridity()` | Aridity index at point locations from the Global Aridity Index GeoTIFF |
| `odisp_fit()` | Maximum likelihood fit; `fixed = list(gamma = 0)` gives a constant range |
| `odisp_profile()` | Profile log-likelihood for one parameter |
| `predict()` | Range, exposure, prevalence and zero-exposure indicator at new locations |
| `odisp_simulate()`, `simulate()` | Simulate survey outcomes, e.g. for checking that parameters are recovered |

## Computation

Exposure and its derivatives are computed in a single compiled pass over the location-node pairs. The pairs are sorted by distance within each location, and the loop for a location stops at the first node beyond the current range, so only pairs inside the range are visited. The gradient is exact.

For the Sudan data (900 villages, about 150,000 nodes and 1.5 million pairs within 100 km), one likelihood evaluation takes about 0.03 s and a complete fit, including three starting values and the Hessian, about 20 s.

`aggregate_nodes()` gives a further reduction. Merging nodes onto 4 km cells, within classes of flow, roughly halves the number of pairs and changes the estimates very little.

## Bounding the dispersal range

Without further constraints, the fitted range tends to grow to whatever `max_range` allows, because a wide kernel can act as a spatial smoother. Setting `phi_max` bounds the range smoothly:

$$
\log\phi = \log\phi_{\max} - \mathrm{softplus}\big(\log\phi_{\max} - \eta_0 - \gamma\log\mathrm{AI}\big).
$$

This equals $\eta_0 + \gamma\log\mathrm{AI}$ well below the bound and approaches $\phi_{\max}$ smoothly. Pairs then only need computing up to `phi_max`, which also makes fitting much faster:

```r
pairs <- river_pairs(villages, nodes, max_range = 40)
fit <- odisp_fit(y ~ 1, ~ log(DIS_AV_CMS) + I(log(DIS_AV_CMS)^2), data = villages,
                 nodes = nodes, pairs = pairs, max_range = 40, phi_max = 40)
```

## Practical notes

- **Choosing `max_range`.** It must exceed the largest fitted dispersal range. Within the fit, a penalty keeps the range below `max_range`, and `odisp_fit()` warns if the range reaches it. When that happens, refit with a larger value. Fits in which the range hits the cap are flagged (`range_capped`) and their log-likelihoods should not be compared.
- **Positive counts far from any river** have zero likelihood without a background term. If such locations exist, either keep smaller streams in the network or set `background = TRUE`.
- **Breeding covariates** should use fixed transformations, for example `log(DIS_AV_CMS)` or `I(log(DIS_AV_CMS)^2)`. Avoid data-dependent bases such as `poly()`, which would differ between the fit and the prediction.
- **Standard errors** assume binomial variation only. Real survey data usually show extra variation between villages, so the standard errors from this version are optimistic. The latent-field version will address this.
