# Pool independent chains for one dataset and mark model; examine mixing for
# every village as well as global parameters. Output spatial coordinates stay
# in survey/node row order. Arguments: inputs.rds chain_directory marks
library(posterior)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("supply inputs, chain directory and mark model")
a <- readRDS(args[1])
out <- args[2]
marks <- match.arg(args[3], c("unit", "two_point"))
seeds <- c(8106, 8206, 8306, 8406)
chains <- lapply(seeds, function(seed) readRDS(file.path(out, sprintf("%s_%d.rds", marks, seed))))
n_keep <- nrow(chains[[1]]$draws)
links <- vapply(chains, function(f) if (is.null(f$settings$link)) "cloglog" else f$settings$link, "")
if (length(unique(links)) != 1L) stop("cannot pool chains with different village links")
for (fit in chains) {
  stopifnot(identical(fit$marks, marks), identical(fit$cells$node, chains[[1]]$cells$node),
            nrow(fit$draws) == n_keep, ncol(fit$prevalence_draws) == nrow(a$villages),
            max(abs(colMeans(fit$prevalence_draws) - fit$fitted$prevalence)) < 1e-12)
}
summarise <- function(x) {
  quantiles <- unname(quantile(x, c(.025, .5, .975)))
  c(mean = mean(x), lo = quantiles[1], median = quantiles[2], hi = quantiles[3],
    rhat = rhat(x), bulk_ess = ess_bulk(x), tail_ess = ess_tail(x), mcse = mcse_mean(x))
}
global <- do.call(rbind, lapply(colnames(chains[[1]]$draws), function(k) {
  x <- vapply(chains, function(f) f$draws[, k], numeric(n_keep))
  data.frame(parameter = k, as.list(summarise(x)), check.names = FALSE)
}))
location <- do.call(rbind, lapply(seq_len(nrow(a$villages)), function(i) {
  x <- vapply(chains, function(f) f$prevalence_draws[, i], numeric(n_keep))
  c(village = i, summarise(x))
}))
cells <- chains[[1]]$cells
for (k in setdiff(names(cells), "node")) {
  cells[[k]] <- rowMeans(vapply(chains, function(f) f$cells[[k]], numeric(nrow(cells))))
}
cells$x_km <- a$nodes$x[cells$node]
cells$y_km <- a$nodes$y[cells$node]
prefix <- file.path(out, marks)
write.csv(global, paste0(prefix, "_parameters.csv"), row.names = FALSE)
predictors <- c(if (length(a$source_covariates)) paste0("theta.",a$source_covariates),
                if (length(a$survey_covariates)) paste0("beta.",a$survey_covariates))
if (length(predictors)) {
  effects <- global[match(predictors, global$parameter), c("parameter","lo","median","hi","rhat","bulk_ess")]
  survey_scale <- if (links[1] == "logit") "village odds ratio" else "village hazard ratio"
  effects$scale <- ifelse(startsWith(effects$parameter,"theta."), "source rate ratio", survey_scale)
  effects[c("lo","median","hi")] <- exp(effects[c("lo","median","hi")])
  write.csv(effects, paste0(prefix, "_covariate_effects.csv"), row.names = FALSE)
}
write.csv(location, paste0(prefix, "_village_risk.csv"), row.names = FALSE)
saveRDS(cells, paste0(prefix, "_cells.rds"))
runs <- data.frame(seed = seeds, elapsed_s = vapply(chains, "[[", 0, "elapsed_seconds"),
                   count_acceptance = vapply(chains, function(f) f$acceptance[["count"]], 0))
if ("relocate" %in% names(chains[[1]]$acceptance))
  runs$relocation_acceptance <- vapply(chains, function(f) f$acceptance[["relocate"]], 0)
write.csv(runs, paste0(prefix, "_runs.csv"), row.names = FALSE)
print(global, digits = 5)
cat("Village median/max R-hat:", median(location[, "rhat"], na.rm = TRUE), max(location[, "rhat"], na.rm = TRUE),
    "; above 1.01:", sum(location[, "rhat"] > 1.01, na.rm = TRUE),
    "; undefined:", sum(is.na(location[, "rhat"])), "; max MCSE:", max(location[, "mcse"], na.rm=TRUE), "\n")
