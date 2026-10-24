###############################################################################
# Joint fit of the tail-up model (range estimated) started from the latent
# mode of its best profile fit. The posterior of the latent effects of the
# tail-up model has several modes, and a fit started from zero ends in a
# lower one (23_directed_profiles.R, step "joint"); here the fit starts from
# the mode of the profile fit at 400 km and, at each evaluation, also from
# the best mode found so far (multistart).
###############################################################################
library(odispers)
library(Matrix)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
inp <- load_result("inputs.rds"); net <- load_result("inputs_full_network.rds")
graph <- river_graph(net$nodes, vertex_spacing = 5)
fits <- readRDS(file.path(out_dir, "directed_tailup.rds"))
b <- fits[[which.max(sapply(fits, `[[`, "logLik"))]]
cat("starting from the profile fit at", b$hyper[["range_u"]], "km, logLik", b$logLik, "\n")
xs <- list(beta = b$coefficients[1:2], c = b$coefficients[3:4], theta = b$coefficients[5], u = b$u, z = b$z)
t0 <- Sys.time()
f <- odisp_fit_zil(y ~ log(AI), ~ log(DIS_AV_CMS), inp$villages, net$nodes, directed_graph(graph, "tailup"),
                   net$pairs, phi = phi_main, range_u = b$hyper[["range_u"]], start = b$par, bfgs = TRUE,
                   estimate_range = TRUE, x_start = xs, checkpoint = file.path(out_dir, "ck_tailup_joint_warm.rds"))
f$time <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
f$setup$graph <- NULL; f$fld$Rch <- NULL
print(f); print(f$hyper)
old <- file.path(out_dir, "directed_tailup_joint.rds")
if (file.exists(old)) file.rename(old, file.path(out_dir, "directed_tailup_joint_coldstart.rds"))
saveRDS(f, old)
