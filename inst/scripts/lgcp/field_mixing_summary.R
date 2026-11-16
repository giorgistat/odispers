# Summarize the four monitored field chains and draw traces/autocorrelations.
# Run from the repository root after field_mixing.R has completed all chains:
# Rscript inst/scripts/lgcp/field_mixing_summary.R output_dir
# Diagnostics use all retained draws; only trace displays are thinned.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("supply the monitored-chain output directory")
library(posterior)
output_dir <- args[1]
chains <- lapply(1:4, function(chain) readRDS(file.path(output_dir, paste0("chain_", chain, ".rds"))))
burnin <- chains[[1]]$burnin
retained <- seq.int(burnin + 1L, nrow(chains[[1]]$field_draws))
variables <- colnames(chains[[1]]$field_draws)
field_array <- array(NA_real_, c(length(retained), length(chains), length(variables)),
                     dimnames = list(NULL, paste0("chain", 1:4), variables))
for (chain in seq_along(chains)) {
  stopifnot(identical(colnames(chains[[chain]]$field_draws), variables),
            identical(chains[[chain]]$vertex_metadata, chains[[1]]$vertex_metadata),
            identical(chains[[chain]]$burnin, burnin))
  field_array[, chain, ] <- chains[[chain]]$field_draws[retained, ]
}

diagnostics <- as.data.frame(summarise_draws(as_draws_array(field_array),
  "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"))
lag_diagnostics <- lapply(seq_along(variables), function(variable) {
  autocorrelations <- vapply(seq_along(chains), function(chain) {
    as.numeric(acf(field_array[, chain, variable], lag.max = 100, plot = FALSE)$acf)[c(2, 11, 101)]
  }, numeric(3))
  data.frame(variable = variables[variable],
             lag1_mean = mean(autocorrelations[1, ]),
             lag10_mean = mean(autocorrelations[2, ]),
             lag100_mean = mean(autocorrelations[3, ]))
})
diagnostics <- merge(diagnostics, do.call(rbind, lag_diagnostics), by = "variable", sort = FALSE)
diagnostics <- diagnostics[match(variables, diagnostics$variable), ]
metadata <- chains[[1]]$vertex_metadata
diagnostics$kind <- "field average or contrast"
diagnostics$kind[diagnostics$variable %in% metadata$variable[metadata$in_footprint]] <- "vertex in footprint"
diagnostics$kind[diagnostics$variable %in% metadata$variable[!metadata$in_footprint]] <- "vertex outside footprint"
write.csv(diagnostics, file.path(output_dir, "field_diagnostics.csv"), row.names = FALSE)
write.csv(metadata, file.path(output_dir, "monitored_vertices.csv"), row.names = FALSE)
write.csv(chains[[1]]$partition_sizes, file.path(output_dir, "partition_sizes.csv"), row.names = FALSE)

# Also retain joint-parameter diagnostics to distinguish field behavior from
# width/count/regression disagreement in these same new runs.
parameter_array <- array(NA_real_, c(length(retained), 4, ncol(chains[[1]]$fit$draws)),
                         dimnames = list(NULL, paste0("chain", 1:4), colnames(chains[[1]]$fit$draws)))
for (chain in 1:4) parameter_array[, chain, ] <- chains[[chain]]$fit$draws
parameter_diagnostics <- as.data.frame(summarise_draws(as_draws_array(parameter_array),
  "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"))
write.csv(parameter_diagnostics, file.path(output_dir, "parameter_diagnostics.csv"), row.names = FALSE)

# The first two plotted functionals were fixed in advance. The local example
# is explicitly selected as the worst R-hat among the monitored footprint sites.
local <- diagnostics[diagnostics$kind == "vertex in footprint", ]
worst_vertex <- local$variable[which.max(local$rhat)]
shown_variables <- c("field_mean_footprint", "field_north_minus_south", worst_vertex)
plot_titles <- c("Mean field in survey footprint", "North minus south mean field",
                 paste("Largest local R-hat:", sub("g_vertex_", "vertex ", worst_vertex)))
displayed <- seq(1L, length(retained), by = 5L)
correlations <- array(NA_real_, c(201, 4, length(shown_variables)),
                      dimnames = list(NULL, NULL, shown_variables))
for (variable in shown_variables) for (chain in 1:4) {
  correlations[, chain, variable] <- as.numeric(acf(field_array[, chain, variable],
                                                   lag.max = 200, plot = FALSE)$acf)
}
plot_data <- list(variables = shown_variables, titles = plot_titles,
                  diagnostics = diagnostics[match(shown_variables, diagnostics$variable), ],
                  traces = field_array[displayed, , shown_variables, drop = FALSE],
                  sweeps = retained[displayed], lags = 0:200,
                  autocorrelations = correlations, retained_draws = length(retained) * 4L)
saveRDS(plot_data, file.path(output_dir, "field_plot_data.rds"))
source("inst/scripts/lgcp/field_mixing_figure.R")
draw_lgcp_field_mixing(plot_data, file.path(output_dir, "field_mixing.png"))

# Compact distributions prevent a few well-mixing coordinates from standing in
# for the whole field. These are diagnostics of 160 monitored vertices, not all.
group_summary <- do.call(rbind, lapply(split(diagnostics, diagnostics$kind), function(group) {
  data.frame(kind = group$kind[1], variables = nrow(group),
    rhat_median = median(group$rhat), rhat_max = max(group$rhat),
    rhat_above_1_01 = sum(group$rhat > 1.01),
    bulk_ess_min = min(group$ess_bulk), bulk_ess_median = median(group$ess_bulk),
    bulk_ess_max = max(group$ess_bulk),
    lag1_median = median(group$lag1_mean), lag100_median = median(group$lag100_mean))
}))
write.csv(group_summary, file.path(output_dir, "field_group_summary.csv"), row.names = FALSE)
shown <- seq(1L, nrow(chains[[1]]$field_draws), by = 5L)
trace_data <- do.call(rbind, lapply(1:4, function(chain) {
  data.frame(chain = chain, sweep = shown, chains[[chain]]$field_draws[shown, , drop = FALSE],
             check.names = FALSE)
}))
write.csv(trace_data, file.path(output_dir, "field_traces_display.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "summary_session.txt"))
print(diagnostics[diagnostics$variable %in% shown_variables, ])
print(diagnostics[diagnostics$kind == "field average or contrast", ])
print(group_summary)
