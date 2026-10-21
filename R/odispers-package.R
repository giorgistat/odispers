#' odispers: river-network dispersal models for onchocerciasis prevalence
#'
#' The main functions are \code{\link{river_nodes}}, which turns a river
#' network into quadrature nodes, \code{\link{river_pairs}}, which finds the
#' nodes within reach of each location, and \code{\link{odisp_fit}}, which
#' fits the model by maximum likelihood.
#'
#' @keywords internal
#' @importFrom stats coef vcov logLik fitted nobs predict simulate
"_PACKAGE"
