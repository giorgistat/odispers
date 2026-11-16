# Independent check of the marked posterior on one 2-km cell shared by two
# villages. Integrate over the field and spread, and sum over total/high sites.
library(odispers)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("supply at most one village link")
link <- if (length(args)) match.arg(args[1], c("cloglog", "logit")) else "cloglog"
village_probability <- function(exposure) {
  if (link == "logit") exposure / (1 + exposure) else -expm1(-exposure)
}
library(sf)
crs <- "+proj=aeqd +lat_0=0 +lon_0=30 +datum=WGS84 +units=km +no_defs"
river <- st_sf(HYRIV_ID = 1L, NEXT_DOWN = 0L, geometry = st_sfc(
  st_linestring(rbind(c(0, 0), c(2, 0))), crs = crs))
nodes <- river_nodes(river, spacing = 3, crs = crs)
graph <- river_graph(nodes)
v <- st_as_sf(data.frame(x = c(1, 1.5), y_coord = 0, n = c(14, 10), y = c(8, 2)),
              coords = c("x", "y_coord"), crs = crs)
K <- wendland(c(0, 0.5)/5) / (5*2/3)
v_field <- as.numeric(graph$A %*% solve(as.matrix(wm_precision(graph, 2, 0.55))) %*%
                        Matrix::t(graph$A))
mu <- function(z) 1.6 * exp(sqrt(v_field) * z)
N <- 0:80
pN <- vapply(N, function(n) integrate(function(z) dnorm(z) * dpois(n, mu(z)),
                  -10, 10, abs.tol = 1e-15)$value, 0)
tail_bound <- integrate(function(z) dnorm(z) * ppois(max(N), mu(z), lower.tail = FALSE),
                        -10, 10, abs.tol = 1e-15)$value + 2*pnorm(-10)
integrate_spread <- function(steps) {
  d <- (seq_len(steps) - 0.5) / steps
  norm <- mean_n <- mean_d <- 0
  for (n in N) for (h in 0:n) {
    mass <- (n-h)*(1-d) + h*(1+d)
    weight <- pN[n+1L] * dbinom(h, n, 0.5) *
      dbinom(v$y[1], v$n[1], village_probability(K[1]*mass)) *
      dbinom(v$y[2], v$n[2], village_probability(K[2]*mass))
    norm <- norm + sum(weight)
    mean_n <- mean_n + n * sum(weight)
    mean_d <- mean_d + sum(d * weight)
  }
  c(total_sites = mean_n/norm, mark_spread = mean_d/norm)
}
ref <- integrate_spread(1024)
stopifnot(max(abs(ref-integrate_spread(512))) < 1e-5, tail_bound < 1e-9)
chains <- lapply(c(901, 902, 903, 904), function(seed) {
  set.seed(seed)
  odisp_fit_lgcp(y~0, ~0, v, nodes, graph, phi=5, range=2, sigma_u=.55,
                 source_rate=.8, marks="two_point", field_update="block_slice", link=link,
                 n_iter=2300L, burnin=300L)$draws
})
for (key in names(ref)) {
  batches <- unlist(lapply(chains, function(x) tapply(x[,key], rep(1:20, each=100),mean)))
  mcse <- sd(batches)/sqrt(length(batches))
  print(c(reference=unname(ref[key]), sampled=mean(batches), mcse=mcse))
  stopifnot(abs(mean(batches)-ref[key]) < 4*mcse)
}
