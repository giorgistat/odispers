# Small REMO pilot. From an installed copy of odispers, run
# Rscript remo_pilot.R <output-directory> reference
# Repeat with half_rate, fine_cells or slice for the sensitivity checks.
library(odispers)
library(sf)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("supply an output directory and a scenario name")
scenario <- match.arg(args[2], c("reference", "half_rate", "fine_cells", "slice"))
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

phi <- 20
range <- 25
sigma_u <- 0.5
spacing <- if (scenario == "fine_cells") 1 else 2
source_rate <- if (scenario == "half_rate") 0.02 else 0.04
seeds <- c(8106, 8206, 8306, 8406)
field_update <- if (scenario == "slice") "slice" else "laplace"
n_iter <- 5000L
burnin <- 1000L
crs <- "+proj=aeqd +lat_0=5 +lon_0=29 +datum=WGS84 +units=km +no_defs"

# The geographic selection uses coordinates only. Include the eight nearest
# surveys to the first Sudan record, including the anchor itself. Retain whole
# reaches within 30 km of any selected village for a 10-km field buffer.
data(remo)
data(sudan_rivers)
v <- st_as_sf(subset(remo, country == "SUDAN"), coords = c("long", "lat"), crs = 4326)
v <- st_transform(v, crs)
xy <- st_coordinates(v)
distance <- sqrt(rowSums(sweep(xy, 2, xy[1, ], "-")^2))
v <- v[head(order(distance), 8), ]
v$y <- round(v$y)
rivers <- st_transform(sudan_rivers, crs)
near <- lengths(st_is_within_distance(rivers, v, dist = phi + 10)) > 0
nodes <- river_nodes(rivers[near, ], spacing = spacing, crs = crs)
graph <- river_graph(nodes, vertex_spacing = 5)
pairs <- river_pairs(v, nodes, phi)
print(graph)
print(pairs)
chains <- lapply(seeds, function(seed) {
  set.seed(seed)
  odisp_fit_lgcp(y ~ 1, ~ 0, v, nodes, graph, pairs, phi = phi, range = range,
                 sigma_u = sigma_u, source_rate = source_rate, prior_sd = 1.5,
                 block_size = 24, n_iter = n_iter, burnin = burnin, field_update = field_update)
})
# Per-chain summaries refer to cells whose midpoints are within kernel range.
# Keep the node coordinates so the retained cell summaries can be located.
saveRDS(list(chains = chains, seeds = seeds, nodes = nodes, graph = graph, surveys = v,
             phi = phi, range = range, sigma_u = sigma_u, source_rate = source_rate),
        file.path(out, paste0(scenario, ".rds")))
diagnostics <- do.call(rbind, lapply(colnames(chains[[1]]$draws), function(name) {
  draws <- vapply(chains, function(f) f$draws[, name], numeric(n_iter - burnin))
  data.frame(parameter = name, mean = mean(draws),
             lo = unname(quantile(draws, 0.025)), hi = unname(quantile(draws, 0.975)),
             split_rhat = posterior::rhat(draws), bulk_ess = posterior::ess_bulk(draws))
}))
accept <- t(vapply(chains, function(f) f$acceptance, c(field = 0, count = 0)))
print(diagnostics, digits = 5)
print(accept)
locations <- data.frame(survey = rownames(v), positives = v$y, examined = v$n,
                        fitted = rowMeans(vapply(chains, function(f) f$fitted$prevalence,
                                                numeric(nrow(v)))))
write.csv(diagnostics, file.path(out, paste0(scenario, "_chains.csv")), row.names = FALSE)
write.csv(locations, file.path(out, paste0(scenario, "_villages.csv")), row.names = FALSE)
write.csv(data.frame(seed = seeds, accept), file.path(out, paste0(scenario, "_acceptance.csv")),
          row.names = FALSE)
