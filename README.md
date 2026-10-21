# odispers

River-network dispersal models for onchocerciasis prevalence.

`odispers` fits exposure models for prevalence survey data in which vector breeding takes place on a river network, exposure at each location accumulates from all river segments within reach, and the distance over which rivers spread infection depends on aridity. Prevalence is exactly zero at locations beyond the reach of any river, so areas without transmission follow from the model rather than being specified in advance.

This first version implements the model **without latent random fields**, fitted by maximum likelihood. The Gaussian field on the river network and the Gaussian process on the plane will follow.

## Model

For a survey location $x$ with $n$ people examined and $y$ found positive,

$$
w(s) = \exp\{z(s)^\top\theta\}, \qquad
E(x) = \int_\Gamma \frac{1}{\phi(x)}\, k\!\left(\frac{\lVert x - s\rVert}{\phi(x)}\right) w(s)\, ds,
$$

$$
\log\phi(x) = \eta_0 + \gamma \log \mathrm{AI}(x), \qquad
Y \sim \mathrm{Bin}\Big(n,\; 1 - \exp\big[-\{E(x) + \epsilon\}\, e^{d(x)^\top\beta}\big]\Big),
$$

where

- $\Gamma$ is the river network and $z(s)$ are characteristics of the river at $s$, such as log discharge;
- $k(t) = (1-t)_+^4(1+4t)$ is the Wendland kernel, which is zero beyond $t = 1$, so a river contributes to $E(x)$ only within the dispersal range $\phi(x)$;
- $\mathrm{AI}(x)$ is the aridity index (precipitation divided by reference evapotranspiration), so $\gamma > 0$ means a shorter range where it is drier;
- $d(x)$ are covariates on the plane, and $\epsilon \ge 0$ is an optional small background exposure.

The integral is computed by quadrature over short pieces of the network.

## Installation

```r
# install.packages("remotes")
remotes::install_github("giorgistat/odispers")
```

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

## Main functions

| Function | Purpose |
|---|---|
| `river_nodes()` | Splits the river network into quadrature nodes with their lengths and attributes |
| `river_pairs()` | Finds the nodes within `max_range` of each location |
| `extract_aridity()` | Aridity index at point locations from the Global Aridity Index GeoTIFF |
| `odisp_fit()` | Maximum likelihood fit; `fixed = list(gamma = 0)` gives a constant range |
| `odisp_profile()` | Profile log-likelihood for one parameter |
| `predict()` | Range, exposure, prevalence and zero-exposure indicator at new locations |
| `odisp_simulate()`, `simulate()` | Simulate survey outcomes, e.g. for checking that parameters are recovered |

## Practical notes

- **Choosing `max_range`.** It must exceed the largest fitted dispersal range. Within the fit, a penalty keeps the range below `max_range`, and `odisp_fit()` warns if the range reaches it. When that happens, refit with a larger value. Fits in which the range hits the cap are flagged (`range_capped`) and their log-likelihoods should not be compared.
- **Positive counts far from any river** have zero likelihood without a background term. If such locations exist, either keep smaller streams in the network or set `background = TRUE`.
- **Breeding covariates** should use fixed transformations, for example `log(DIS_AV_CMS)` or `I(log(DIS_AV_CMS)^2)`. Avoid data-dependent bases such as `poly()`, which would differ between the fit and the prediction.
- **Standard errors** assume binomial variation only. Real survey data usually show extra variation between villages, so the standard errors from this version are optimistic. The latent-field version will address this.
