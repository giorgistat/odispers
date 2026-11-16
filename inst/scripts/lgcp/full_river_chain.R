# One chain on the full Sudan survey. Input contains villages, nodes, graph
# and 40-km pairs in the same order. See the LGCP README for preparation.
# Arguments: inputs.rds output_directory unit|two_point seed [simulated_y.rds] [--relocate|--extended] [--intercept-only]
library(odispers)
args <- commandArgs(trailingOnly = TRUE)
intercept_only <- "--intercept-only" %in% args
grid_run <- "--grid" %in% args
link_arg <- grep("^--link=", args, value=TRUE)
if (length(link_arg)>1L) stop("supply a single village link")
link <- if (length(link_arg)) match.arg(sub("^--link=", "", link_arg), c("cloglog","logit")) else "cloglog"
range_arg <- grep("^--field-range=",args,value=TRUE)
if (length(range_arg)>1L) stop("supply a single --field-range in km")
field_range <- if (length(range_arg)) as.numeric(sub("^--field-range=","",range_arg)) else 25
if (!is.finite(field_range) || field_range<=0) stop("invalid field range")
extended <- "--extended" %in% args
long_run <- grid_run || extended || "--relocate" %in% args
args <- args[!args %in% c("--relocate","--extended","--intercept-only","--grid",range_arg,link_arg)]
if (!length(args) %in% 4:5) stop("supply inputs, output directory, marks, seed and optional simulated counts")
a <- readRDS(args[1])
out <- args[2]
marks <- match.arg(args[3], c("unit", "two_point"))
seed <- as.integer(args[4])
if (length(args) == 5L) {
  simulated <- readRDS(args[5])
  if (length(simulated) != nrow(a$villages)) stop("one simulated outcome per survey is required")
  a$villages$y <- simulated
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!intercept_only && (!length(a$survey_covariates) || !length(a$source_covariates)))
  stop("prepare environmental predictors using full_river_covariates.R, or specify --intercept-only for a baseline")
survey_formula <- if (!intercept_only) reformulate(a$survey_covariates, response = "y") else y ~ 1
source_formula <- if (!intercept_only) reformulate(a$source_covariates, intercept = FALSE) else ~ 0
set.seed(seed)
elapsed <- system.time(fit <- odisp_fit_lgcp(survey_formula, source_formula, a$villages, a$nodes, a$graph, a$pairs,
                         phi = if (is.null(a$phi)) 40 else a$phi, range = field_range, sigma_u = 0.5, source_rate = 0.04,
                         prior_sd = 1.5, block_size = if (long_run) 1L else 24L,
                         n_iter = if (grid_run) 6000L else if (extended) 14000L else if (long_run) 8000L else 5000L,
                         burnin = if (grid_run) 1000L else if (extended) 2500L else if (long_run) 1500L else 1000L, relocate = long_run,
                         field_update = "block_slice", marks = marks, keep_prevalence = TRUE, link = link))
fit$elapsed_seconds <- unname(elapsed["elapsed"])
fit$seed <- seed
fit$formula <- survey_formula
fit$breeding <- source_formula
fit$survey_input <- if (length(args) == 5L) args[5] else args[1]
destination <- file.path(out, sprintf("%s_%d.rds", marks, seed))
tmp <- tempfile(tmpdir = out)
saveRDS(fit, tmp, compress = FALSE)
if (!file.rename(tmp, destination)) stop("could not publish chain output")
print(fit)
print(elapsed)
