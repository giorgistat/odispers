# Joint LGCP fit using all training surveys. Test outcomes have zero likelihood weight.
library(odispers)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==5L)
out <- args[1]
model <- match.arg(args[2],c("baseline","river_climate","occurrence","both","flexible_occurrence","survey_design","unit_occurrence","exposure_occurrence"))
fold <- as.integer(args[3])
seed <- as.integer(args[4])
sweeps <- as.integer(args[5])
a <- readRDS(file.path(out,"inputs.rds"))
v <- a$villages
# Sensitivity conditional on the recorded survey design. Keep the original stratum
# when masking the held-out response and trial count below.
v$examined30 <- as.integer(v$n==30L)
v$examined50 <- as.integer(v$n==50L)
if (fold>0L) {
  held <- a$spatial_fold == fold
  v$y[held] <- v$n[held] <- 0L
}
river_terms <- a$source_covariates
if (model %in% c("river_climate","both")) river_terms <- c(river_terms,"river_dry")
arguments <- list(formula=y~z_logAI,breeding=reformulate(river_terms,intercept=FALSE),
                  data=v,nodes=a$nodes,graph=a$graph,pairs=a$pairs,phi=40,range=25,
                  sigma_u=.5,source_rate=.04,prior_sd=1.5,marks="two_point",link="logit",
                  block_size=1L,relocate=TRUE,field_update="block_slice",keep_prevalence=TRUE,
                  n_iter=sweeps,burnin=floor(sweeps*.2))
if (model %in% c("occurrence","both","unit_occurrence","exposure_occurrence")) arguments$occurrence <- ~climate
if (model=="exposure_occurrence") arguments$occurrence_exposure <- TRUE
if (model=="unit_occurrence") arguments$marks <- "unit"
if (model=="flexible_occurrence") arguments$occurrence <- ~climate + I(pmin(climate + .75, 0))
if (model=="survey_design") arguments$occurrence <- ~climate + I(pmin(climate + .75, 0)) + examined30 + examined50
set.seed(seed)
elapsed <- system.time(fit <- do.call(odisp_fit_lgcp, arguments))
# Keep the fit compact without retaining the embedded data passed by do.call().
fit$call <- NULL
fit$seed <- seed
fit$fold <- fold
fit$training <- which(v$n>0)
fit$elapsed_seconds <- unname(elapsed["elapsed"])
dir <- file.path(out,model,paste0("fold",fold))
dir.create(dir,showWarnings=FALSE,recursive=TRUE)
temp <- tempfile(tmpdir=dir)
saveRDS(fit,temp,compress=FALSE)
stopifnot(file.rename(temp,file.path(dir,paste0("chain_",seed,".rds"))))
print(fit)
cat("Elapsed seconds",fit$elapsed_seconds,"\n")
