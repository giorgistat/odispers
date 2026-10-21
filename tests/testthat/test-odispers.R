test_that("Wendland kernel has the right values and support", {
  expect_equal(wendland(0), 1)
  expect_equal(wendland(c(1, 1.5)), c(0, 0))
  expect_equal(wendland(0.5), 0.5^4 * 3)
})

test_that("quadrature preserves river length", {
  s <- toy_setup(m = 20)
  len_true <- sum(as.numeric(sf::st_length(sf::st_transform(s$rivers, attr(s$nodes, "crs")))))
  expect_equal(sum(s$nodes$len), len_true, tolerance = 1e-6)
  expect_true(all(s$nodes$len <= 1 + 1e-9))
})

test_that("exposure on the bank of a long straight river does not depend on phi", {
  river <- sf::st_sf(geometry = sf::st_sfc(sf::st_linestring(rbind(c(29, 0), c(31, 0))), crs = 4326))
  pt <- sf::st_sf(y = 0L, n = 1L, AI = 1, geometry = sf::st_sfc(sf::st_point(c(30, 0)), crs = 4326))
  nodes <- river_nodes(river, spacing = 0.1)
  pairs <- river_pairs(pt, nodes, max_range = 100)
  des <- odisp_design(y ~ 1, ~ 0, pt, nodes, pairs)
  E <- sapply(c(10, 20, 40), function(phi)
    odispers:::odisp_eval(c("beta.(Intercept)" = 0, eta0 = log(phi), gamma = 0), des, FALSE)$E)
  ck <- integrate(function(v) wendland(abs(v)), -1, 1)$value
  expect_equal(E, rep(ck, 3), tolerance = 1e-3)
})

test_that("analytic gradient matches numerical gradient", {
  s <- toy_setup(m = 150)
  s$villages$y <- rbinom(150, 40, 0.2)
  des <- odisp_design(y ~ 1, ~ log(flow), s$villages, s$nodes, s$pairs)
  for (bg in c(FALSE, TRUE)) {
    obj <- odispers:::odisp_objective(des, bg, list())
    x <- c(-1, 0.3, log(30), 0.5, if (bg) log(0.01))
    num <- vapply(seq_along(x), function(k) {
      h <- 1e-6; e <- replace(numeric(length(x)), k, h)
      (obj$fn(x + e) - obj$fn(x - e)) / (2 * h)
    }, numeric(1))
    expect_equal(obj$gr(x), num, tolerance = 1e-4)
  }
})

test_that("parameters are recovered from simulated data", {
  s <- toy_setup(m = 600, seed = 2)
  s$villages$y <- 0L
  des <- odisp_design(y ~ 1, ~ log(flow), s$villages, s$nodes, s$pairs)
  truth <- c("beta.(Intercept)" = 0.5, "theta.log(flow)" = 0.3, eta0 = log(25), gamma = 0.6)
  set.seed(3)
  s$villages$y <- odisp_simulate(des, truth)[, 1]
  fit <- odisp_fit(y ~ 1, ~ log(flow), data = s$villages, nodes = s$nodes, pairs = s$pairs)
  est <- coef(fit)[names(truth)]
  se <- sqrt(diag(vcov(fit)))[names(truth)]
  expect_true(all(abs(est - truth) < 3.5 * se))
  expect_s3_class(summary(fit), "summary.odisp_fit")
})

test_that("prevalence is exactly zero beyond reach and predict works", {
  s <- toy_setup(m = 400, seed = 4)
  s$villages$y <- 0L
  des <- odisp_design(y ~ 1, ~ 0, s$villages, s$nodes, s$pairs)
  truth <- c("beta.(Intercept)" = 0, eta0 = log(15), gamma = 0)
  set.seed(5)
  s$villages$y <- odisp_simulate(des, truth)[, 1]
  fit <- odisp_fit(y ~ 1, ~ 0, data = s$villages, nodes = s$nodes, pairs = s$pairs,
                   fixed = list(gamma = 0))
  expect_equal(unname(coef(fit)["gamma"]), 0)
  far <- sf::st_sf(AI = 1, geometry = sf::st_sfc(sf::st_point(c(31.5, 2.0)), crs = 4326))
  pr <- predict(fit, far, s$nodes)
  expect_true(pr$zero_exposure)
  expect_equal(pr$prevalence, 0)
})
