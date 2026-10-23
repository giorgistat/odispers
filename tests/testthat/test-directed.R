toy_graph <- function(q = c(9, 1, 10)) {
  # two tributaries (ids 1, 2) flowing into reach 3, an outlet; discharge q
  nodes <- data.frame(line = rep(1:3, c(3, 2, 4)), x = 0, y = 0,
                      len = rep(c(2, 1.5, 3), c(3, 2, 4)),
                      HYRIV_ID = rep(1:3, c(3, 2, 4)), NEXT_DOWN = rep(c(3, 3, 0), c(3, 2, 4)),
                      DIS_AV_CMS = rep(q, c(3, 2, 4)))
  class(nodes) <- c("odisp_network", "data.frame")
  river_graph(nodes, vertex_spacing = 5)
}

test_that("tail-down field is stationary with exponential covariance in river distance", {
  g <- directed_graph(toy_graph(), "taildown")
  S <- as.matrix(solve(wm_precision(g, range = 10, sigma = 1.5)))
  expect_equal(diag(S), rep(1.5^2, g$nv), tolerance = 1e-8)
  k <- 2 / 10
  # source of reach 1 (vertex 1) and source of reach 2 (vertex 2): river distance 6 + 3
  expect_equal(S[1, 2], 1.5^2 * exp(-k * 9), tolerance = 1e-8)
  # source of reach 1 and the outlet: 6 + 12
  expect_equal(S[1, g$nv], 1.5^2 * exp(-k * 18), tolerance = 1e-8)
})

test_that("tail-up field: tributaries independent, flow-weighted mixing below the confluence", {
  g <- directed_graph(toy_graph(q = c(9, 1, 10)), "tailup")
  S <- as.matrix(solve(wm_precision(g, range = 10)))
  k <- 2 / 10
  ns <- g$n_src; e <- function(i) ns + i      # latent value at the downstream end of edge i
  # edges: reach 1 -> 1, 2; reach 2 -> 3; reach 3 -> 4, 5, 6
  expect_equal(diag(S)[seq_len(ns + 3)], rep(1, ns + 3), tolerance = 1e-8)   # above the confluence
  expect_equal(S[e(2), e(3)], 0, tolerance = 1e-12)                          # tributaries independent
  # start of reach 3 is 0.9 * end of reach 1 + 0.1 * end of reach 2
  r4 <- exp(-k * 4)
  expect_equal(S[e(2), e(4)], 0.9 * r4, tolerance = 1e-8)
  expect_equal(S[e(3), e(4)], 0.1 * r4, tolerance = 1e-8)
  expect_equal(S[e(4), e(4)], r4^2 * (0.9^2 + 0.1^2) + 1 - r4^2, tolerance = 1e-8)
  # the projection interpolates between the two ends of each edge
  expect_equal(as.vector(Matrix::rowSums(g$A)), rep(1, 9))
})

test_that("sum of tail-up and tail-down fields splits the variance", {
  g <- directed_graph(toy_graph(), "tailup+taildown")
  Q <- wm_precision(g, range = 10, mix = 0.3)
  S <- as.matrix(Matrix::solve(Q))
  expect_equal(dim(Q), c(g$nv, g$nv))
  expect_equal(S[1, 1], 0.3, tolerance = 1e-8)                   # tail-up source
  expect_equal(S[g$nv, g$nv], 0.7, tolerance = 1e-8)             # tail-down vertex
  expect_equal(S[1, g$nv], 0, tolerance = 1e-12)                 # independent components
})
