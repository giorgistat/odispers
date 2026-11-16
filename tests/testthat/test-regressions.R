# Small deterministic network: two connected reaches and eight nearby villages.
regression_network <- function() {
  crs <- "+proj=aeqd +lat_0=0 +lon_0=30 +datum=WGS84 +units=km +no_defs"
  rivers <- sf::st_sf(HYRIV_ID = 1:2, NEXT_DOWN = c(2, 0), flow = c(1, 4),
                     habitat = factor(c("A", "B")),
                     geometry = sf::st_sfc(
                       sf::st_linestring(rbind(c(0, 0), c(10, 0))),
                       sf::st_linestring(rbind(c(10, 0), c(20, 0))), crs = crs))
  nodes <- river_nodes(rivers, spacing = 5, crs = crs)
  villages <- sf::st_as_sf(data.frame(x = c(2, 4, 6, 8, 12, 14, 16, 18), h = 2,
                                     y = c(0, 3, 1, 5, 0, 8, 2, 10), n = 40, AI = 0.5,
                                     group = factor(rep(c("A", "B"), 4))),
                          coords = c("x", "h"), crs = crs)
  list(v = villages, nodes = nodes, graph = river_graph(nodes, vertex_spacing = 5),
       pairs = river_pairs(villages, nodes, 30), phi = 15)
}

test_that("missing and infinite breeding covariates are rejected before exposure evaluation", {
  s <- regression_network()
  for (bad in c(NA_real_, NaN, Inf, 0)) {
    nodes <- s$nodes; nodes$flow[2] <- bad
    expect_error(odisp_design(y ~ 1, ~ log(flow), s$v, nodes, s$pairs), "missing values|finite")
  }
  d <- odisp_design(y ~ 1, ~ log(flow), s$v, s$nodes, s$pairs)
  expect_equal(nrow(d$Z), length(d$len))
  # A malformed pre-built design must also fail safely at the compiled boundary.
  d$Z <- d$Z[-2, , drop = FALSE]
  expect_error(odispers:::odisp_eval(c(0, 0, log(15), 0), d, FALSE), "one row per node")
})

test_that("prediction ignores observed responses and uses no trial-count placeholders", {
  s <- regression_network()
  f <- odisp_fit(I(y) ~ group, ~ log(flow), s$v, s$nodes, s$pairs,
                 fixed = list(eta0 = log(s$phi), gamma = 0))
  pr <- predict(f, s$v, s$nodes)
  expect_equal(pr$prevalence, fitted(f), tolerance = 1e-12)
  newdata <- s$v; newdata$y <- newdata$n <- NULL
  expect_equal(predict(f, newdata, s$nodes), pr)
  newdata$y <- NA_real_; newdata$n <- NA_real_
  expect_equal(predict(f, newdata, s$nodes), pr)
})

test_that("prediction covers new ranges and rejects incomplete supplied pairs", {
  s <- regression_network(); crs <- attr(s$nodes, "crs")
  nodes <- structure(data.frame(x = c(0, 40), y = 0, len = 1),
                     class = c("odisp_network", "data.frame"), crs = crs)
  training <- sf::st_sf(y = 10, n = 40, AI = 0.5,
                        geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = crs))
  f <- odisp_fit(y ~ 1, ~ 0, training, nodes, max_range = 30,
                 fixed = list(eta0 = log(40), gamma = 1))
  newdata <- sf::st_sf(AI = 2, geometry = sf::st_sfc(sf::st_point(c(80, 0)), crs = crs))
  pr <- predict(f, newdata, nodes)
  # Only the node 40 km away contributes to the 80 km kernel.
  expected_E <- (1 - 0.5)^4 * (1 + 4 * 0.5) / (2 / 3 * 80)
  expect_equal(pr$phi, 80)
  expect_equal(pr$E, expected_E, tolerance = 1e-12)
  expect_gt(pr$prevalence, 0)
  expect_false(pr$zero_exposure)
  expect_equal(predict(f, newdata, nodes, river_pairs(newdata, nodes, 80)), pr)
  expect_error(predict(f, newdata, nodes, river_pairs(newdata, nodes, 30)), "too small")
  f$coefficients["eta0"] <- log(30); newdata$AI <- 1
  expect_no_error(predict(f, newdata, nodes, river_pairs(newdata, nodes, 30)))
})

test_that("prediction preserves prevalence and breeding factor encodings", {
  s <- regression_network()
  contrasts(s$v$group) <- stats::contr.sum(2)
  f <- odisp_fit(y ~ group, ~ habitat, s$v, s$nodes, s$pairs,
                 fixed = list(eta0 = log(s$phi), gamma = 0))
  pr <- predict(f, s$v, s$nodes)
  newdata <- s$v; newdata$group <- factor(as.character(newdata$group), levels = c("B", "A"))
  nodes <- s$nodes; nodes$habitat <- factor(nodes$habitat, levels = c("B", "A"))
  expect_equal(predict(f, newdata, nodes), pr)
  single <- newdata[2, ]; single$group <- factor("B")
  expect_equal(predict(f, single, nodes)$prevalence, pr$prevalence[2])
  newdata$group <- as.character(newdata$group)
  expect_equal(predict(f, newdata, nodes), pr)
  newdata$group[1] <- "new level"
  expect_error(predict(f, newdata, nodes), "new level")
  far <- single; sf::st_geometry(far) <- sf::st_sfc(sf::st_point(c(200, 0)), crs = sf::st_crs(far))
  expect_equal(predict(f, far, nodes)$E, 0)
})

test_that("latent model predictions retain factor encoding and validate breeding inputs", {
  s <- regression_network()
  fits <- list(
    odisp_fit_u(y ~ group, ~ habitat, s$v, s$nodes, s$graph, s$pairs,
                 phi = s$phi, field = "none", trace = FALSE),
    odisp_fit_zi(y ~ group, ~ habitat, s$v, s$nodes, s$graph, s$pairs,
                  phi = s$phi, range_u = 50, trace = FALSE,
                  fixed = list(sigma_u = 0.6, sigma_z = 0.4, alpha = 1, gamma = 0.4, delta = 0)),
    odisp_fit_zil(y ~ group, ~ habitat, s$v, s$nodes, s$graph, s$pairs,
                   phi = s$phi, range_u = 50, trace = FALSE,
                   fixed = list(sigma_u = 0.6, sigma_z = 0.4, a_E = 0.8, gamma = 0.4)))
  newdata <- s$v; newdata$group <- factor(newdata$group, levels = c("B", "A"))
  nodes <- s$nodes; nodes$habitat <- factor(nodes$habitat, levels = c("B", "A"))
  bad_nodes <- nodes; bad_nodes$habitat[1] <- NA
  for (f in fits) {
    pr <- predict(f, s$v, s$nodes, s$graph)
    expect_equal(predict(f, newdata, nodes, s$graph), pr)
    single <- newdata[2, ]; single$group <- factor("B")
    expect_equal(predict(f, single, nodes, s$graph)$prevalence, pr$prevalence[2])
    expect_error(predict(f, newdata, bad_nodes, s$graph), "missing values")
    expect_error(predict(f, newdata, nodes, s$graph, river_pairs(newdata, nodes, 5)), "too small")
    far <- single; sf::st_geometry(far) <- sf::st_sfc(sf::st_point(c(200, 0)), crs = sf::st_crs(far))
    expect_equal(predict(f, far, nodes, s$graph)$E, 0)
  }
})

test_that("saved model terms preserve prediction bases and old numeric fits remain usable", {
  s <- regression_network(); s$v$z <- seq_len(nrow(s$v))
  f <- odisp_fit(y ~ poly(z, 2), ~ log(flow), s$v, s$nodes, s$pairs,
                 fixed = list(eta0 = log(s$phi), gamma = 0))
  pr <- predict(f, s$v, s$nodes)
  expect_equal(predict(f, s$v[2, ], s$nodes)$prevalence, pr$prevalence[2])
  f <- odisp_fit(y ~ z, ~ log(flow), s$v, s$nodes, s$pairs,
                 fixed = list(eta0 = log(s$phi), gamma = 0))
  pr <- predict(f, s$v, s$nodes)
  f$design$prevalence_spec <- f$design$breeding_spec <- NULL
  expect_equal(predict(f, s$v, s$nodes), pr)
  f <- odisp_fit(y ~ group, ~ log(flow), s$v, s$nodes, s$pairs,
                 fixed = list(eta0 = log(s$phi), gamma = 0))
  f$design$prevalence_spec <- NULL
  expect_error(predict(f, s$v, s$nodes), "refit")
})

test_that("constant breeding works in latent and zero-inflated models", {
  s <- regression_network()
  for (breeding in list(~ 0, ~ 1)) {
    fits <- list(
      odisp_fit_u(y ~ 1, breeding, s$v, s$nodes, s$graph, s$pairs,
                   phi = s$phi, field = "none", trace = FALSE),
      odisp_fit_zi(y ~ 1, breeding, s$v, s$nodes, s$graph, s$pairs,
                    phi = s$phi, range_u = 50, trace = FALSE,
                    fixed = list(sigma_u = 0.6, sigma_z = 0.4, alpha = 1, gamma = 0.4, delta = 0)),
      odisp_fit_zil(y ~ 1, breeding, s$v, s$nodes, s$graph, s$pairs,
                     phi = s$phi, range_u = 50, trace = FALSE,
                     fixed = list(sigma_u = 0.6, sigma_z = 0.4, a_E = 0.8, gamma = 0.4)))
    for (f in fits) {
      expect_true(is.finite(f$logLik))
      expect_false(any(grepl("^theta", names(f$coefficients))))
      expect_equal(dim(f$vcov), rep(length(f$coefficients), 2))
      expect_true(all(is.finite(predict(f, s$v, s$nodes, s$graph)$prevalence)))
    }
  }
})

test_that("summed directed fields honour distinct fixed ranges", {
  s <- regression_network(); graph <- directed_graph(s$graph, "tailup+taildown")
  f <- odisp_fit_zil(y ~ 1, ~ log(flow), s$v, s$nodes, graph, s$pairs,
                     phi = s$phi, range_u = 50, range_up = 200, trace = FALSE,
                     fixed = list(sigma_u = 0.6, sigma_z = 0.4, a_E = 0.8, gamma = 0.4, mix = 0.3))
  expect_equal(unname(f$hyper[c("range_u", "range_up")]), c(50, 200))
  expect_equal(f$fld$R, wm_precision(graph, 50, mix = 0.3, range_up = 200))
  expect_true(is.finite(f$logLik))
})

test_that("integrated fitting rejects truncated distance pairs before constructing exposure", {
  s <- regression_network()
  short <- river_pairs(s$v, s$nodes, 5)
  expect_error(odispers:::odisp_u_setup(y ~ 1, ~ log(flow), s$v, s$nodes,
    s$graph, short, s$phi, "n", "AI", 10), "max_range.*too small")
  supplied <- odispers:::odisp_u_setup(y ~ 1, ~ log(flow), s$v, s$nodes,
    s$graph, s$pairs, s$phi, "n", "AI", 10)
  automatic <- odispers:::odisp_u_setup(y ~ 1, ~ log(flow), s$v, s$nodes,
    s$graph, NULL, s$phi, "n", "AI", 10)
  expect_equal(supplied$K, automatic$K)
  expect_equal(supplied$used, automatic$used)
})

test_that("the integrated fitter retains exact parameters when returning a cached best mode", {
  s <- regression_network()
  evaluated <- NULL
  # Return a previously cached inferior point to exercise the best-mode fallback.
  testthat::local_mocked_bindings(odisp_newton_max = function(f, x, ...) {
    f(x)
    candidate <- x + c(.123456789123, -.012345678912, .023456789123, -.001234567891)
    values <- c(f(candidate), f(x))
    evaluated <<- if (values[1] > values[2]) candidate else x
    if (values[1] > values[2]) x else candidate
  }, .package = "odispers")
  fit <- odisp_fit_zil(y ~ 1, ~ log(flow), s$v, s$nodes, s$graph, s$pairs,
    phi = s$phi, range_u = 50, start = c(log(.6), log(.4), .8, .4), trace = FALSE)
  expect_identical(fit$par, evaluated)
  expect_identical(unname(fit$hyper[c("sigma_u", "sigma_z")]), exp(unname(fit$par[1:2])))
  expect_identical(unname(fit$hyper[c("a_E", "gamma")]), unname(fit$par[3:4]))
  state <- list(beta = fit$coefficients[1], c = fit$coefficients[2],
    theta = fit$coefficients[3], u = fit$u, z = fit$z)
  ev <- odispers:::odisp_zil_eval(state, fit$setup, fit$fld)
  expect_equal(exp(ev$a), fit$E, tolerance = 1e-12)
  expect_equal(plogis(ev$d), fit$pi, tolerance = 1e-12)
})

test_that("HPC cross-validation aligns original folds to retained villages", {
  # Load helper definitions only, without running the data preparation or fits.
  script <- system.file("scripts", "hpc", "taildown_covariates.R", package = "odispers")
  helpers <- new.env(parent = asNamespace("odispers"))
  wanted <- c("subset_setup", "gh", "g2", "g2x1", "g2x2", "g2w", "pmf_dE", "cv_zil")
  for (expr in parse(script)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) eval(expr, helpers)
  }
  s <- regression_network()
  far <- s$v[1, ]; far$y <- 0
  sf::st_geometry(far) <- sf::st_sfc(sf::st_point(c(200, 0)), crs = sf::st_crs(far))
  v <- rbind(s$v[1:3, ], far, s$v[4:8, ])
  f <- odisp_fit_zil(y ~ 1, ~ log(flow), v, s$nodes, s$graph, phi = s$phi,
                     range_u = 50, trace = FALSE,
                     fixed = list(sigma_u = 0.6, sigma_z = 0.4, a_E = 0.8, gamma = 0.4))
  folds <- rep(1:2, length.out = nrow(v))
  # This label disappears entirely when the unreachable village is omitted.
  folds[4] <- 3L
  pr <- helpers$cv_zil(f, folds)
  expect_equal(vapply(pr, `[[`, integer(1), "village"), f$setup$has)
  expect_equal(vapply(pr, `[[`, numeric(1), "y"), f$setup$y)
  expect_equal(vapply(pr, function(x) sum(x$pmf), numeric(1)), rep(1, length(pr)), tolerance = 1e-10)
  # Reindexing an equivalent setup must leave its predictive distributions unchanged.
  compact <- f; compact$setup$m_all <- length(f$setup$has)
  compact$setup$has <- seq_along(f$setup$has)
  reference <- helpers$cv_zil(compact, folds[f$setup$has])
  expect_equal(lapply(pr, `[[`, "pmf"), lapply(reference, `[[`, "pmf"), tolerance = 1e-10)
  expect_error(helpers$cv_zil(f, folds[-1]), "original village")
})
