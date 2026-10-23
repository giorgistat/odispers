###############################################################################
# Profile of the approximate log-likelihood over the share of variance phi of
# the tail-up part in the tail-up + tail-down field (21_fit_directed.R), with
# the other parameters at the estimates of the tail-down model (phi = 0).
# phi = 0 is the tail-down model, phi = 1 the tail-up model.
###############################################################################
library(odispers)
library(Matrix)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
inp <- load_result("inputs.rds"); net <- load_result("inputs_full_network.rds")
graph <- river_graph(net$nodes, vertex_spacing = 5)
td <- readRDS(file.path(out_dir, "directed_taildown.rds"))[["800"]]
dg <- directed_graph(graph, "tailup+taildown")
S <- odispers:::odisp_u_setup(y ~ log(AI), ~ log(DIS_AV_CMS), inp$villages, net$nodes, dg, net$pairs,
                              phi_main, "n", "AI", 10)
S$m_c <- td$setup$m_c
h <- td$hyper; nu <- dg$nv_up
phis <- c(0.01, 0.05, 0.1, 0.2, 0.35, 0.5)
x <- list(beta = td$coefficients[1:2], c = td$coefficients[3:4], theta = td$coefficients[5],
          u = c(rep(0, nu), td$u), z = td$z)
out <- data.frame(phi = phis, logLik = NA)
for (k in seq_along(phis)) {
  R <- wm_precision(dg, 800, 1, mix = phis[k])
  fld <- list(field = "u", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE),
              sig_u = h[["sigma_u"]], sig_z = h[["sigma_z"]], aE = h[["a_E"]], gamma = h[["gamma"]])
  r <- odispers:::odisp_zil_laplace(fld, S, x); x <- r$x
  out$logLik[k] <- r$logLik
  cat(sprintf("phi %.2f  logLik %.3f\n", phis[k], r$logLik))
}
out <- rbind(data.frame(phi = 0, logLik = td$logLik), out)
print(out)
saveRDS(out, file.path(out_dir, "profile_phi.rds"))
