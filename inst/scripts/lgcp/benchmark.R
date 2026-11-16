# Time a short full-network chain with identical settings and random seed.
# This measures computation, not convergence or effective samples per second.
# Run from the repository root:
# Rscript inst/scripts/lgcp/benchmark.R inputs.rds library_dir output_prefix
# Use separate R processes for the old and new installed native libraries.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("supply inputs.rds, an installed R library and an output prefix")
library(odispers, lib.loc = args[2])
inputs <- readRDS(args[1])
output_prefix <- args[3]

seed <- 111L
iterations <- 300L
burnin <- 100L
repetitions <- 4L
timings <- numeric(repetitions)
for (repetition in seq_len(repetitions)) {
  set.seed(seed)
  gc()
  elapsed <- system.time(fit <- odisp_fit_lgcp(
    y ~ z_logAI, breeding = reformulate(inputs$source_covariates, intercept = FALSE),
    data = inputs$villages, nodes = inputs$nodes, graph = inputs$graph,
    pairs = inputs$pairs, phi = 60, range = 400, sigma_u = .5,
    source_rate = .04, prior_sd = 1.5, marks = "two_point", link = "logit",
    occurrence = ~ climate, occurrence_exposure = TRUE,
    block_size = 1L, relocate = TRUE, field_update = "block_slice",
    keep_prevalence = TRUE, n_iter = iterations, burnin = burnin,
    kernel = "gaussian", width = 10, estimate_width = TRUE,
    width_bounds = c(1, 240)))
  timings[repetition] <- unname(elapsed["elapsed"])
  cat("Repetition", repetition, ":", timings[repetition], "seconds\n")
  flush.console()
}

fit$call <- NULL
saveRDS(fit, paste0(output_prefix, "_fit.rds"))
write.csv(data.frame(repetition = seq_len(repetitions), seconds = timings,
                     iterations = iterations, burnin = burnin, seed = seed),
          paste0(output_prefix, "_timings.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = paste0(output_prefix, "_session.txt"))
