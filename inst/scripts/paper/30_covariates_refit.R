###############################################################################
# Refits of the tail-down model with RiverATLAS breeding covariates selected
# by 29_covariates_explore.R, with sigma_u, sigma_z, a_E and gamma
# re-estimated (the range is fixed at the joint tail-down estimate, about
# 1017 km, to save time). Covariates are standardised over the nodes that
# enter the likelihood.
# Usage: Rscript 30_covariates_refit.R <slope|inundation|forest|slope_inund_forest>
###############################################################################
library(odispers)
library(Matrix)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
model <- commandArgs(trailingOnly = TRUE)[1]
inp <- load_result("inputs.rds"); v <- inp$villages
net <- load_result("inputs_full_network.rds"); nodes <- net$nodes
data(sudan_river_covariates)
rc <- as.data.frame(sudan_river_covariates)[match(nodes$HYRIV_ID, sudan_river_covariates$HYRIV_ID), ]
td <- readRDS(file.path(out_dir, "directed_taildown_joint.rds"))
used <- td$setup$used
std <- function(x) (x - mean(x[used])) / sd(x[used])
nodes$slope <- std(log1p(rc$sgr_dk_rav))
nodes$inundation <- std(log1p(rc$inu_pc_ult))
nodes$forest <- std(rc$for_pc_use)
covs <- list(slope = "slope", inundation = "inundation", forest = "forest",
             slope_inund_forest = c("slope", "inundation", "forest"))[[model]]
breeding <- stats::as.formula(paste("~ log(DIS_AV_CMS) +", paste(covs, collapse = " + ")))
g <- directed_graph(river_graph(nodes, vertex_spacing = 5), "taildown")
cf <- td$coefficients
x0 <- list(beta = cf[1:2], c = cf[3:4], theta = c(cf[5], rep(0, length(covs))), u = td$u, z = td$z)
t0 <- Sys.time()
f <- odisp_fit_zil(y ~ log(AI), breeding, v, nodes, g, net$pairs, phi = phi_main,
                   range_u = td$hyper[["range_u"]], start = td$par[1:4], bfgs = TRUE, x_start = x0,
                   multistart = FALSE, checkpoint = file.path(out_dir, sprintf("ck_cov_%s.rds", model)))
f$time <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
f$setup$graph <- NULL; f$fld$Rch <- NULL
print(f)
cat(sprintf("change in logLik against the tail-down model: %.2f\n", f$logLik - td$logLik))
saveRDS(f, file.path(out_dir, sprintf("cov_refit_%s.rds", model)))
