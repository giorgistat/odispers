###############################################################################
# Full model with both spatial random effects and the village effect:
#   w(s) = exp{theta log Q(s) + u(s)},
#   log T(x) = beta0 + beta1 log AI(x) + S(x) + Z(x),
# u on the river network (in the exposure) and S on the plane (in the
# intensity of transmission), range of dispersal 40 km. The ranges of u and
# S are fitted on a small grid around the best ranges of the models with one
# field (10_fit_s_vs_u.R). Compared with the models with one field.
###############################################################################

library(odispers)
library(sf)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

inp <- load_result("inputs.rds")
v <- inp$villages; nodes <- inp$nodes
pairs <- river_pairs(v, nodes, phi_main)
graph <- river_graph(nodes, vertex_spacing = 5)

# pairs of ranges (km): u along the network, S in straight-line distance
range_grid <- rbind(c(800, 400), c(800, 800), c(800, 1600), c(400, 800))
fits <- odisp_fit_u(y ~ log(AI), ~ log(DIS_AV_CMS), v, nodes, graph, pairs, phi = phi_main,
                    field = "uS", range = range_grid, z = TRUE)
prof <- attr(fits, "profile")
cat("\nModel with u, S and Z:\n"); print(prof, row.names = FALSE)
best <- fits[[which.max(prof$logLik)]]
print(best)

# comparison with the models with one field, if 10_fit_s_vs_u.R has been run
f10 <- file.path(out_dir, "s_vs_u_fits.rds")
if (file.exists(f10)) {
  sm <- readRDS(f10)$summary
  cat("\nBest model of each kind (Laplace log-likelihood):\n")
  print(rbind(sm[, c("model", "description", "logLik")],
              data.frame(model = "ET", description = "u in E + S in T + village effect Z",
                         logLik = round(best$logLik, 2))), row.names = FALSE)
}
saveRDS(list(fits = fits, profile = prof, best = best), file.path(out_dir, "u_s_z_fits.rds"))
