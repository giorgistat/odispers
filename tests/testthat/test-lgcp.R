test_that("cell lengths enter Poisson mass once and shared events each contribute one kernel", {
  a <- lgcp_test_data()
  S <- lgcp_setup(y ~ 0, ~ 0, a$data, a$nodes, a$graph, NULL, 5, 2, 0.55,
                  0.8, "n", "line", 2)
  kernel <- c(1, 0.9^4 * 1.4, 0) / (5 * 2/3)
  expect_equal(exp(S$log_mass), 1.6)
  expect_equal(as.numeric(S$K %*% 2), 2 * kernel)
  # For two shared surveys the probability of zero cases is one joint
  # Poisson transform. Separate village transforms give a different answer.
  expected <- exp(1.6 * expm1(-sum(S$n * kernel)))
  ll <- vapply(0:26, function(k) lgcp_survey(k * kernel, rep(0, 3), 0*S$y, S$n)$loglik, 0)
  expect_equal(sum(stats::dpois(0:26, 1.6) * exp(ll)), expected, tolerance = 1e-12)
  expect_gt(expected, prod(exp(1.6 * expm1(-S$n * kernel))))
  expect_equal(lgcp_survey(1, -1000, 1, 1)$loglik, -1000)
})

test_that("conditional Laplace uses the Poisson score and observed curvature", {
  a <- lgcp_test_data()
  S <- lgcp_setup(y ~ 0, ~ 1, a$data, a$nodes, a$graph, NULL, 5, 2, 0.55,
                  0.8, "n", "line", 2)
  counts <- 3L
  lap <- lgcp_field_laplace(counts, S)
  x <- lap$mode
  step <- 1e-4
  ev <- function(x) lgcp_field_log_density(x, counts, S)
  for (k in seq_along(x)) {
    e <- diag(length(x))[, k] * step
    expect_lt(abs((ev(x + e) - ev(x - e)) / (2 * step)), 2e-6)
    curvature <- (2 * ev(x) - ev(x + e) - ev(x - e)) / step^2
    expect_equal(as.numeric(lap$H[k, k]), curvature, tolerance = 1e-6)
  }
})

test_that("retained draws reconstruct both villages and preserve structural zeros", {
  withr::local_seed(1106)
  a <- lgcp_test_data()
  f <- odisp_fit_lgcp(y ~ 0, ~ 0, a$data, a$nodes, a$graph, phi = 5,
                      range = 2, sigma_u = 0.55, source_rate = 0.8,
                      n_iter = 100, burnin = 10, keep_counts = TRUE)
  expect_type(f$count_draws, "integer")
  expect_equal(unname(f$draws[, "total_sites"]), as.numeric(f$count_draws))
  expect_true(all(f$count_draws >= 1L))
  expect_equal(f$cells$mean_count, mean(f$count_draws))
  k <- c(1, 0.9^4 * 1.4, 0) / (5 * 2/3)
  expect_equal(f$fitted$exposure, mean(f$count_draws) * k)
  expect_equal(f$fitted$prevalence,
               colMeans(-expm1(-outer(as.numeric(f$count_draws), k))))
  expect_s3_class(f, "odisp_lgcpfit")
  expect_equal(fitted(f), f$fitted$prevalence)
})

test_that("the sampler uses offsets and both coefficient formulas", {
  withr::local_seed(191)
  a <- lgcp_test_data()
  a$nodes$z <- 0.2
  a$data$age <- c(0, 0.1, -0.2)
  S <- lgcp_setup(y ~ age + offset(age), ~ z + offset(z), a$data, a$nodes,
                  a$graph, NULL, 5, 2, 0.55, 0.8, "n", "area", 2)
  expect_equal(S$survey_offset, a$data$age)
  expect_equal(exp(S$log_mass), 1.6 * exp(0.2))
  expect_equal(as.numeric(S$K[1, 1]), 7 / (pi * 5^2))
  f <- odisp_fit_lgcp(y ~ age, ~ z, a$data, a$nodes, a$graph, phi = 5,
                      n_iter = 35, burnin = 20)
  expect_named(as.data.frame(f$draws), c("beta.(Intercept)", "beta.age",
                                         "theta.(Intercept)", "theta.z", "total_sites"))
  expect_true(all(is.finite(f$draws)))
})

test_that("slice updates use the joint counts and a draw from the permuted graph prior", {
  withr::local_seed(921)
  a <- lgcp_test_data()
  f <- odisp_fit_lgcp(y ~ 0, ~ 1, a$data, a$nodes, a$graph, phi = 5,
                      n_iter = 75, burnin = 15, field_update = "slice", keep_counts = TRUE)
  expect_true(is.na(f$acceptance[["field"]]))
  expect_true(all(f$count_draws > 0))
  expect_equal(f$fitted$exposure[1], mean(f$count_draws) / (5 * 2 / 3))

  # A star puts the high-degree first vertex late in the sparse elimination.
  Q <- Matrix::sparseMatrix(i = rep(1, 4), j = 2:5, x = 0.4, symmetric = TRUE) +
    Matrix::Diagonal(5, 2)
  ch <- Matrix::Cholesky(Q, LDL = FALSE)
  z <- withr::with_seed(359, stats::rnorm(5))
  x <- withr::with_seed(359, lgcp_gaussian_draw(ch, 5))
  expect_equal(sum(x * as.numeric(Q %*% x)), sum(z^2), tolerance = 1e-12)
})

test_that("an incomplete kernel search or positive unreachable count is rejected", {
  a <- lgcp_test_data()
  short <- river_pairs(a$data, a$nodes, max_range = 2)
  expect_error(odisp_fit_lgcp(y ~ 1, data = a$data, nodes = a$nodes, graph = a$graph,
                              pairs = short, phi = 5), "max_range")
  a$data$y[3] <- 1
  expect_error(odisp_fit_lgcp(y ~ 1, data = a$data, nodes = a$nodes, graph = a$graph,
                              phi = 5), "beyond the reach")
})

test_that("the local count updates agree with full recomputation for unit and marked sites", {
  S <- list(K = Matrix::Matrix(rbind(c(0.8, 0, 0), c(0.2, 0.5, 0), c(0, 1e-18, 1)),
                               sparse = TRUE), y = c(1, 2, 1), n = c(5, 4, 3))
  S$kernel_tree <- lgcp_tree_layout(S$K)
  mu <- c(0.8, 2.1, 1.2)
  predictor <- c(-1, .5, 0)
  # Independent reference: compute every village likelihood for every trial.
  reference <- function(initial, types, spread) {
    marked <- length(types) > 0
    ll <- function(n, h) {
      mass <- if (marked) (n - h) * (1-spread) + h * (1+spread) else n
      E <- as.numeric(S$K %*% mass)
      odds <- E * exp(predictor)
      probability <- if (S$link == "logit") odds/(1 + odds) else -expm1(-odds)
      value <- sum(stats::dbinom(S$y, S$n, probability, log = TRUE))
      if (!is.null(S$activity_sign)) {
        g <- plogis(S$gate_base + S$gate_slope*E)
        value <- value + log(g[1]) + log1p(-g[2])
      }
      value
    }
    counts <- initial
    order <- sample.int(3)
    proposals <- accepts <- 0L
    for (ids in list(order[1:2], order[3])) {
      trial <- counts
      high <- types
      for (j in ids) {
        trial[j] <- stats::rpois(1, mu[j])
        if (marked) high[j] <- stats::rbinom(1, trial[j], .5)
      }
      if (identical(counts, trial) && identical(high, types)) next
      proposals <- proposals + 1L
      if (log(stats::runif(1)) < ll(trial, high) - ll(counts, types)) {
        counts <- trial
        types <- high
        accepts <- accepts + 1L
      }
    }
    list(counts = counts, high = types, proposals = proposals, accepts = accepts)
  }
  for (link in c("cloglog", "logit")) for (marked in c(FALSE, TRUE)) for (coupled in c(FALSE, TRUE)) {
    S$link <- link
    S$activity_sign <- if (coupled) c(1L,-1L,0L)
    S$gate_base <- if (coupled) c(-.7,1,.3)
    S$gate_slope <- if (coupled) -2.5
    S$y <- if (coupled) c(1,0,0) else c(1,2,1)
    S$n <- if (coupled) c(5,0,0) else c(5,4,3)
    counts <- c(1L, 1L, 0L)
    high <- if (marked) c(0L, 1L, 0L) else integer()
    d <- if (marked) 0.7 else 0
    for (seed in 1:25) {
      slow <- withr::with_seed(seed, reference(counts, high, d))
      fast <- withr::with_seed(seed, lgcp_count_sweep(counts, mu, predictor, S, 2L, high, d))
      expect_identical(fast[names(slow)], slow)
      mass <- if (marked) (fast$counts-fast$high)*(1-d) + fast$high*(1+d) else fast$counts
      expect_equal(fast$exposure, as.numeric(S$K %*% mass))
      counts <- fast$counts
      high <- fast$high
    }
  }
})

test_that("the marked draws reconstruct village prevalence and cell strength", {
  a <- lgcp_test_data()
  withr::local_seed(91106)
  f <- odisp_fit_lgcp(y ~ 0, ~ 0, a$data, a$nodes, a$graph, phi = 5, range = 2,
                      sigma_u = 0.55, source_rate = 0.8, n_iter = 100, burnin = 10,
                      field_update = "block_slice", marks = "two_point", keep_counts = TRUE,
                      keep_prevalence = TRUE)
  d <- f$draws[, "mark_spread"]
  H <- as.numeric(f$high_draws)
  N <- as.numeric(f$count_draws)
  expect_true(all(H >= 0 & H <= N & d >= 0 & d <= 1))
  mass <- (1-d) * (N-H) + (1+d) * H
  k <- c(1, .9^4*1.4, 0)/(5*2/3)
  expect_equal(f$prevalence_draws, -expm1(-outer(mass, k)))
  expect_equal(f$cells$mean_marked_mass, mean(mass))
  expect_equal(f$cells$mean_high_count, mean(H))
  expect_equal(fitted(f), colMeans(f$prevalence_draws))
})

test_that("graph blocks respect precision and projected cell boundaries", {
  # A branching graph with an isolated vertex and several bands.
  from <- c(1:31, 13:20)
  to <- c(2:32, 33:40)
  for (shift in c(0, 6)) {
    b <- lgcp_graph_blocks(from, to, 41, 12, shift)
    same <- b$colour[from] == b$colour[to]
    expect_true(any(!same))
    expect_equal(b$group[from[same]], b$group[to[same]])
    expect_true(all(b$group[from[!same]] != b$group[to[!same]]))
  }
})

test_that("site relocation targets the marked shared-village posterior on two cells", {
  withr::local_seed(6243)
  K <- matrix(c(1, .2, .3, .7), 2, 2)
  mu <- c(.6, .9)
  d <- .7
  y <- c(2, 0); n <- c(3, 2)
  S <- list(K = Matrix::Matrix(K, sparse = TRUE), y = y, n = n)
  S$kernel_tree <- lgcp_tree_layout(S$K)
  # Independently split each Poisson process into two equally likely types.
  # With per-type mean <= .45, the omitted mass beyond 10 of a type is <1e-11.
  states <- expand.grid(low1 = 0:10, high1 = 0:10, low2 = 0:10, high2 = 0:10)
  N <- cbind(states$low1 + states$high1, states$low2 + states$high2)
  H <- cbind(states$high1, states$high2)
  E <- ((N-H)*(1-d) + H*(1+d)) %*% t(K)
  w <- stats::dbinom(y[1], n[1], -expm1(-E[,1])) *
    stats::dbinom(y[2], n[2], -expm1(-E[,2]))
  for (j in seq_len(4L)) w <- w * stats::dpois(states[[j]], rep(mu/2, each = 2)[j])
  expect_lt(4 * stats::ppois(10, max(mu)/2, lower.tail = FALSE), 2e-11)
  # Start each independent replicate at a draw from the enumerated posterior.
  # A valid sweep preserves both the count and mark marginal at each cell.
  ids <- sample.int(nrow(states), 1800, replace = TRUE, prob = w)
  differences <- t(vapply(ids, function(i) {
    ct <- as.integer(N[i,]); h <- as.integer(H[i,])
    moved <- lgcp_count_sweep(ct, mu, c(0,0), S, 1L, h, d, c(1L,1L))
    c(moved$counts-ct, moved$high-h)
  }, numeric(4)))
  mcse <- apply(differences, 2, stats::sd)/sqrt(length(ids))
  expect_true(all(abs(colMeans(differences)) < 4.5 * mcse))
})

test_that("blocked field updates support breeding coefficients and the correct Gaussian conditional", {
  withr::local_seed(2419)
  a <- lgcp_test_data()
  a$graph <- river_graph(a$nodes, vertex_spacing = 0.5)
  S <- lgcp_setup(y ~ 0, ~ 1, a$data, a$nodes, a$graph, NULL, 5, 2, .55, .8, "n", "line", 2)
  blocks <- lgcp_block_setup(S, a$graph, shift = 0L, width = 1L)
  C <- solve(as.matrix(wm_precision(a$graph, range = 2, sigma = .55)))
  for (b in blocks) {
    # Dense Gaussian conditional mean/covariance are independent references.
    regression <- C[b$idx,b$other,drop=FALSE] %*% solve(C[b$other,b$other,drop=FALSE])
    expected_cov <- C[b$idx,b$idx,drop=FALSE] - regression %*% C[b$other,b$idx,drop=FALSE]
    expect_equal(as.matrix(Matrix::solve(b$ch, -b$cross)), regression, tolerance = 1e-10)
    expect_equal(as.matrix(Matrix::solve(b$ch, Matrix::Diagonal(length(b$idx)))),
                 expected_cov, tolerance = 1e-10)
  }
  x <- lgcp_block_step(numeric(ncol(S$B)), 1L, S, blocks)
  expect_length(x, a$graph$nv + 1L)
  expect_true(all(is.finite(x)) && x[1] != 0)
})

test_that("the breeding-coefficient update matches a Poisson regression quadrature", {
  withr::local_seed(84107)
  Z <- matrix(c(-1.3, .4, 1.5), ncol=1)
  offset <- log(c(.9, 2.1, 1.4))
  counts <- c(1L, 4L, 3L)
  prior_sd <- 1.5
  density <- function(theta) vapply(theta, function(t) {
    prod(stats::dpois(counts, exp(offset + Z[,1]*t))) * stats::dnorm(t, sd=prior_sd)
  }, 0)
  norm <- stats::integrate(density, -10, 10)$value
  ref <- stats::integrate(function(x) x*density(x), -10, 10)$value / norm
  x <- 0
  draws <- numeric(2500)
  for (i in seq_along(draws)) {
    x <- lgcp_breeding_step(x, Z, offset, counts, prior_sd)
    draws[i] <- x
  }
  batches <- tapply(draws[501:2500], rep(1:40, each=50), mean)
  expect_lt(abs(mean(batches)-ref), 4.5*stats::sd(batches)/sqrt(length(batches)))
})

test_that("tail-down field blocks match the tree covariance at short and long ranges", {
  # A branching tree with projection at edge midpoints. Distances below are
  # an independent covariance reference, including the unconnected tributaries.
  graph <- structure(list(from = c(1L, 2L, 3L, 4L), to = c(3L, 3L, 4L, 5L),
    length = c(2, 3, 4, 1), nv = 5L, type = "wm", node_edge = 1:4,
    node_fr = rep(.5, 4), edge_q = c(1, 1, 2, 2),
    A = Matrix::sparseMatrix(i = rep(1:4, 2), j = c(1, 2, 3, 4, 3, 3, 4, 5),
      x = .5, dims = c(4, 5))), class = "odisp_graph")
  graph <- directed_graph(graph, "taildown")
  distance <- matrix(c(0,5,2,6,7, 5,0,3,7,8, 2,3,0,4,5, 6,7,4,0,1, 7,8,5,1,0), 5)
  for (range in c(2, 800)) {
    Q <- wm_precision(graph, range, .5)
    C <- .5^2 * exp(-2 * distance / range)
    expect_equal(as.matrix(solve(Q)), C, tolerance = 1e-10)
    S <- list(Q = Q, used = 1:4, theta_names = character())
    blocks <- lgcp_block_setup(S, graph, shift = 0L, width = 1L)
    for (b in blocks) {
      regression <- C[b$idx, b$other, drop = FALSE] %*% solve(C[b$other, b$other, drop = FALSE])
      conditional <- C[b$idx, b$idx, drop = FALSE] - regression %*% C[b$other, b$idx, drop = FALSE]
      expect_equal(as.matrix(Matrix::solve(b$ch, -b$cross)), regression, tolerance = 1e-10)
      expect_equal(as.matrix(Matrix::solve(b$ch, Matrix::Diagonal(length(b$idx)))), conditional, tolerance = 1e-10)
    }
  }
  withr::local_seed(8713)
  a <- lgcp_test_data()
  fit <- odisp_fit_lgcp(y ~ 1, ~ 0, a$data, a$nodes, directed_graph(a$graph, "taildown"),
    phi = 5, range = 800, n_iter = 40, burnin = 10, field_update = "block_slice",
    occurrence = ~ 1, occurrence_exposure = TRUE, marks = "two_point", link = "logit")
  expect_true(all(is.finite(fit$draws)))
  expect_equal(fit$fitted$prevalence[3], 0)
})

test_that("cloglog regression curvature agrees with numerical derivatives", {
  h <- c(0, 1e-9, .03, .7, 3, 30)
  y <- c(0, 1, 2, 1, 4, 8); n <- c(5, 5, 5, 5, 5, 8)
  d <- lgcp_cloglog_score(h, y, n)
  for (i in seq_along(h)) {
    ll <- function(eta) lgcp_survey(h[i], eta, y[i], n[i])$loglik
    score <- function(eta) lgcp_cloglog_score(h[i]*exp(eta),y[i],n[i])$score
    step <- 1e-5
    expect_equal(d$score[i], (ll(step)-ll(-step))/(2*step), tolerance=1e-7)
    expect_equal(d$curvature[i], (score(-step)-score(step))/(2*step), tolerance=1e-7)
  }
})

test_that("the local Gaussian update and slice sampling agree with a binomial posterior", {
  withr::local_seed(8751)
  S <- list(X=matrix(c(-1,.2,1.5),ncol=1), survey_offset=c(-.2,0,.3),
            n=c(15,10,18), y=c(1,2,9), prior_sd=1.5)
  E <- c(.4,.15,.3)
  for (link in c("cloglog", "logit")) {
    S$link <- link
    p <- function(x) {
      odds <- E * exp(S$survey_offset + as.numeric(S$X)*x)
      if (link == "logit") odds / (1 + odds) else -expm1(-odds)
    }
    density <- function(x) vapply(x,function(v) prod(stats::dbinom(S$y,S$n,p(v))) *
                                  stats::dnorm(v,sd=S$prior_sd),0)
    ref <- stats::integrate(function(x) x*density(x),-10,10,abs.tol=1e-13)$value /
      stats::integrate(density,-10,10,abs.tol=1e-13)$value
    b <- 0; draws <- numeric(2500)
    for(i in seq_along(draws)) {
      b <- if(i %% 4L) lgcp_beta_newton(b,E,S) else lgcp_beta_step(b,E,S)
      draws[i] <- b
    }
    batches <- tapply(draws[501:2500], rep(1:40,each=50),mean)
    expect_lt(abs(mean(batches)-ref),4.5*stats::sd(batches)/sqrt(length(batches)))
  }
})

test_that("exposure trees retain tiny surviving weights and exact structural zeros", {
  withr::local_seed(8040)
  K <- Matrix::sparseMatrix(i=c(1,1,2),j=c(1,2,1),x=c(1,1e-24,0.7),dims=c(3,2))
  S <- list(K=K,kernel_tree=lgcp_tree_layout(K),n=c(0,0,0),y=c(0,0,0))
  a <- lgcp_count_sweep(c(1L,1L),c(0,20),numeric(3),S,1)
  expect_equal(a$counts[1],0L)
  expect_gt(a$counts[2],0)
  expect_identical(a$exposure, c(a$counts[2]*1e-24,0,0))
  b <- lgcp_count_sweep(a$counts,c(0,0),numeric(3),S,1)
  expect_identical(b$exposure,numeric(3))
})


test_that("logit survey likelihood handles zero exposure and both numerical tails", {
  expect_equal(lgcp_survey(c(0, 2), c(1000, 0), c(0, 1), c(5, 4), "logit")$probability,
               c(0, 2/3))
  expect_equal(lgcp_survey(0, 0, 1, 2, "logit")$loglik, -Inf)
  expect_equal(lgcp_survey(1, -1000, 1, 3, "logit")$loglik, -1000)
  expect_equal(lgcp_survey(1, 1000, 2, 3, "logit")$loglik, -1000)
})

test_that("logit village draws reconstruct from shared site strengths and survey offsets", {
  withr::local_seed(8050)
  a <- lgcp_test_data()
  a$data$shift <- c(-1, .5, 1)
  for (method in c("laplace", "slice", "block_slice")) {
    fit <- odisp_fit_lgcp(y ~ 0 + offset(shift), ~0, a$data, a$nodes, a$graph,
                          phi=5, range=2, sigma_u=.55, source_rate=.8, n_iter=60,
                          burnin=10, marks="two_point", keep_counts=TRUE, keep_prevalence=TRUE,
                          field_update=method, relocate=TRUE, link="logit")
    spread <- fit$draws[, "mark_spread"]
    high <- as.numeric(fit$high_draws)
    strength <- (as.numeric(fit$count_draws) - high)*(1-spread) + high*(1+spread)
    kernel <- c(1, .9^4*1.4, 0)/(5*2/3)
    odds <- outer(strength, kernel * exp(a$data$shift))
    expect_equal(fit$settings$link, "logit")
    expect_equal(fit$prevalence_draws, odds/(1+odds))
  }
})

test_that("one shared occurrence state updates with the full binomial probability of zero", {
  g <- c(.2, .7, .8, .4)
  E <- c(.5, 0, 1, 2)
  n <- c(12, 50, 0, 14)
  y <- c(0, 0, 0, 1)
  for (link in c("cloglog", "logit")) {
    p <- if (link=="logit") E/(1+E) else -expm1(-E)
    zero <- (1-g)+g*(1-p)^n
    exact <- g*(1-p)^n/zero
    exact[y>0] <- 1
    expect_equal(lgcp_activity_probability(qlogis(g),E,rep(0,4),y,n,link),exact)
    expect_gt(zero[1],(1-g[1]*p[1])^n[1])
  }
})

test_that("fixed occurrence offsets enforce inactive survey zeros and preserve positive villages", {
  a <- lgcp_test_data()
  a$data$y[2] <- 0
  a$data$offset_gate <- c(0,-1000,1000)
  withr::local_seed(67100)
  f <- odisp_fit_lgcp(y~0,~0,a$data,a$nodes,a$graph,phi=5,range=2,sigma_u=.55,
                      source_rate=.8,link="logit",occurrence=~0+offset(offset_gate),
                      field_update="block_slice",marks="two_point",n_iter=60,burnin=10,
                      keep_prevalence=TRUE)
  expect_equal(f$fitted$pr_active,c(1,0,1))
  expect_true(all(f$prevalence_draws[,1]>0))
  expect_true(all(f$prevalence_draws[,2:3]==0))
  expect_error(odisp_fit_lgcp(y~1,~0,a$data,a$nodes,a$graph,phi=5,occurrence=y~1),"one-sided")
})
