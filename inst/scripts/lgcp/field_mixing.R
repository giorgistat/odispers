# Monitor the latent Gaussian field in a full-data LGCP chain.
# Run from the repository root with an installed copy of the current package:
# Rscript inst/scripts/lgcp/field_mixing.R inputs.rds library_dir output_dir chain
# Chains 1--4 use the report's seeds, width starts and 8,000/2,000 schedule.
# The local wrapper records field values without changing proposals or RNG calls.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("supply inputs.rds, library_dir, output_dir and chain (1--4)")
library(odispers, lib.loc = args[2])
library(Matrix)

chain <- as.integer(args[4])
stopifnot(length(chain) == 1L, chain %in% 1:4)
seeds <- c(111L, 222L, 333L, 444L)
initial_widths <- c(10, 20, 40, 80)
iterations <- 8000L
burnin <- 2000L
range_km <- 800
cutoff_km <- 60
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
inputs <- readRDS(args[1])
inputs$graph <- directed_graph(inputs$graph, "taildown")
graph <- inputs$graph
used <- sort(unique(inputs$pairs$j[inputs$pairs$r < cutoff_km]))
projection <- graph$A[used, , drop = FALSE]
in_footprint <- which(Matrix::colSums(abs(projection)) > 0)
outside <- setdiff(seq_len(graph$nv), in_footprint)

# Fixed ID quantiles avoid choosing locations after looking at their traces.
spaced_ids <- function(ids, number) {
  ids[unique(round(seq(1, length(ids), length.out = min(number, length(ids)))))]
}
vertices <- c(spaced_ids(in_footprint, 128L), spaced_ids(outside, 32L))
vertex_names <- paste0("g_vertex_", vertices)
vertex_metadata <- data.frame(variable = vertex_names, vertex = vertices,
                              in_footprint = vertices %in% in_footprint)

# Cell-length-weighted field averages probe broad spatial motion. Multiplying
# their weights by the interpolation matrix avoids reconstructing every cell.
cell_length <- inputs$nodes$len
latitude <- inputs$nodes$z_lat[used]
latitude_cut <- median(latitude)
north <- used[latitude >= latitude_cut]
south <- used[latitude < latitude_cut]
field_weights <- function(cells) {
  as.numeric(Matrix::crossprod(graph$A[cells, , drop = FALSE],
                              cell_length[cells] / sum(cell_length[cells])))
}
average_weights <- rbind(
  field_mean_all = field_weights(seq_along(cell_length)),
  field_mean_footprint = field_weights(used),
  field_mean_north = field_weights(north),
  field_mean_south = field_weights(south))
average_weights <- rbind(average_weights,
  field_north_minus_south = average_weights["field_mean_north", ] -
                            average_weights["field_mean_south", ])
monitor_names <- c(rownames(average_weights), vertex_names)

# Each field call completes both colours. The enclosing fitter updates beta
# afterwards; that does not change the recorded field, so this is an end-of-sweep
# field trajectory. A private function environment leaves the package untouched.
package_namespace <- asNamespace("odispers")
ordinary_step <- get("lgcp_block_step", envir = package_namespace)
monitor_environment <- new.env(parent = package_namespace)
field_draws <- matrix(NA_real_, iterations, length(monitor_names),
                     dimnames = list(NULL, monitor_names))
sweep <- 0L
monitor_environment$lgcp_block_step <- function(state, counts, S, blocks) {
  updated <- ordinary_step(state, counts, S, blocks)
  field <- updated[length(S$theta_names) + seq_len(graph$nv)]
  sweep <<- sweep + 1L
  field_draws[sweep, ] <<- c(as.numeric(average_weights %*% field), field[vertices])
  if (sweep %% 500L == 0L) {
    cat("Chain", chain, "sweep", sweep, "at", format(Sys.time()), "\n")
    flush.console()
  }
  updated
}
monitored_fit <- get("odisp_fit_lgcp", envir = package_namespace)
environment(monitored_fit) <- monitor_environment
fit_arguments <- list(
  formula = y ~ z_logAI,
  breeding = reformulate(inputs$source_covariates, intercept = FALSE),
  data = inputs$villages, nodes = inputs$nodes, graph = graph, pairs = inputs$pairs,
  phi = cutoff_km, range = range_km, sigma_u = .5, source_rate = .04,
  prior_sd = 1.5, marks = "two_point", link = "logit",
  occurrence = ~ climate, occurrence_exposure = TRUE,
  block_size = 1L, relocate = TRUE, field_update = "block_slice",
  keep_prevalence = FALSE, n_iter = iterations, burnin = burnin,
  kernel = "gaussian", width = initial_widths[chain], estimate_width = TRUE,
  width_bounds = c(1, 240))

# A short identical-seed replay verifies that recording leaves every returned
# quantity and the final RNG state unchanged before the diagnostic run begins.
short_arguments <- fit_arguments
short_arguments$n_iter <- 12L
short_arguments$burnin <- 4L
set.seed(seeds[chain])
ordinary <- do.call(odisp_fit_lgcp, short_arguments)
ordinary_rng <- .Random.seed
set.seed(seeds[chain])
monitored <- do.call(monitored_fit, short_arguments)
monitored_rng <- .Random.seed
ordinary$call <- monitored$call <- NULL
stopifnot(identical(ordinary, monitored), identical(ordinary_rng, monitored_rng),
          sweep == short_arguments$n_iter, all(is.finite(field_draws[seq_len(sweep), ])))
cat("Monitoring replay: identical fit and RNG state\n")
sweep <- 0L
field_draws[,] <- NA_real_

# The graph partition width is 12 edge steps, shifted by six every other sweep.
# Record its actual connected-component sizes rather than calling it pointwise.
partition_sizes <- do.call(rbind, lapply(c(0L, 6L), function(shift) {
  labels <- get("lgcp_graph_blocks", envir = package_namespace)(
    graph$from, graph$to, graph$nv, 12L, shift)
  size <- tabulate(labels$group)
  data.frame(shift = shift, group = seq_along(size), vertices = size,
             touches_footprint = seq_along(size) %in% labels$group[in_footprint])
}))

set.seed(seeds[chain])
cat("START tail-down 800 km, chain", chain, "initial width", initial_widths[chain], "\n")
elapsed <- system.time(fit <- do.call(monitored_fit, fit_arguments))[["elapsed"]]
stopifnot(sweep == iterations, all(is.finite(field_draws)))
fit$call <- NULL
saveRDS(list(field_draws = field_draws, fit = fit, burnin = burnin,
             seed = seeds[chain], initial_width = initial_widths[chain],
             elapsed_seconds = elapsed, vertex_metadata = vertex_metadata,
             latitude_cut = latitude_cut, partition_sizes = partition_sizes),
        file.path(output_dir, paste0("chain_", chain, ".rds")))
capture.output(sessionInfo(), file = file.path(output_dir, paste0("session_", chain, ".txt")))
cat("DONE chain", chain, "seconds", elapsed, "\n")
