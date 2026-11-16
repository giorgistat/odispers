# Assemble the fixed LGCP design: river-cell masses, sparse village kernels,
# field interpolation and Gaussian precision. Cell order is shared throughout.
lgcp_setup <- function(formula, breeding, data, nodes, graph, pairs, phi, range,
                       sigma_u, source_rate, trials, scaling, prior_sd,
                       kernel = "wendland", width = phi/3) {
  if (!inherits(data, "sf") || !inherits(nodes, "odisp_network") ||
      !inherits(graph, "odisp_graph") || nrow(graph$A) != nrow(nodes))
    stop("supply survey locations, river nodes and a graph built from the same nodes")
  pairs <- odisp_prediction_pairs(data, nodes, pairs, phi)
  inside <- pairs$r < phi
  i <- pairs$i[inside]
  j <- pairs$j[inside]
  used <- sort(unique(j))
  survey <- odisp_model_matrix(formula, sf::st_drop_geometry(data))
  y <- stats::model.response(survey$frame)
  n <- data[[trials]]
  if (!is.numeric(y) || !is.numeric(n) || !is.null(dim(y)) || !is.null(dim(n)) ||
      length(y) != nrow(data) || length(n) != nrow(data) ||
      any(!is.finite(c(y, n))) || any(y != floor(y) | n != floor(n) | y < 0 | n < y))
    stop("survey counts must be integers satisfying 0 <= y <= trials")
  if (any(y > 0 & tabulate(i, nrow(data)) == 0L))
    stop("villages with positives are beyond the reach of every river cell")
  if (!length(used)) stop("no river cells within the survey kernel radius")
  lengths <- nodes$len[used]
  if (any(!is.finite(lengths) | lengths <= 0)) stop("river cell lengths must be positive")
  if (length(breeding) != 2L) stop("breeding must be a one-sided formula")
  habitat <- odisp_model_matrix(breeding, as.data.frame(nodes)[used, , drop = FALSE])
  prior <- wm_precision(graph, range, sigma = sigma_u)
  B <- graph$A[used, , drop = FALSE]
  n_theta <- ncol(habitat$matrix)
  if (n_theta) {
    B <- cbind(Matrix::Matrix(habitat$matrix, sparse = TRUE), B)
    prior <- Matrix::bdiag(Matrix::Diagonal(n_theta, 1 / prior_sd^2), prior)
  }
  scale <- if (scaling == "line") 2 / 3 * phi else pi / 7 * phi^2
  kernel_distance <- NULL
  if (kernel == "wendland") {
    K <- Matrix::sparseMatrix(i = i, j = match(j, used),
                              x = wendland(pairs$r[inside] / phi) / scale,
                              dims = c(length(y), length(used)))
  } else {
    # CSC order is by source cell, then village. Keep explicit zero entries
    # so very narrow widths cannot change the sum-tree's structural layout.
    K <- Matrix::sparseMatrix(i=i, j=match(j,used), x=rep(1,length(i)),
                              dims=c(length(y),length(used)))
    if (length(K@x) != length(i)) stop("Gaussian kernel requires unique village-cell distance pairs")
    kernel_distance <- pairs$r[inside][order(j,i)]
    K@x <- lgcp_gaussian_kernel(kernel_distance,width,phi,scaling)
  }
  if (any(!is.finite(K@x)))
    stop("kernel weights must be finite; check the kernel width and scaling")
  offset <- function(frame) {
    v <- stats::model.offset(frame)
    if (is.null(v)) return(numeric(nrow(frame)))
    if (any(!is.finite(v))) stop("formula offsets must be finite")
    as.numeric(v)
  }
  log_mass <- log(lengths) + log(source_rate) + offset(habitat$frame)
  baseline_mass <- exp(log_mass)
  if (any(!is.finite(baseline_mass) | baseline_mass == 0))
    stop("baseline Poisson cell means must be representable and positive; check source_rate and offsets")
  list(K = K, kernel_tree = lgcp_tree_layout(K), kernel_distance=kernel_distance,
       kernel_cutoff=phi, kernel_scaling=scaling, B = B, Q = Matrix::forceSymmetric(prior), X = survey$matrix,
       breeding_matrix = habitat$matrix, field_projection = graph$A[used, , drop = FALSE],
       used = used, length = lengths, log_mass = log_mass,
       survey_offset = offset(survey$frame), theta_names = colnames(habitat$matrix),
       y = y, n = n, prior_sd = prior_sd)
}
