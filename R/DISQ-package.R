#' Distributed inference for spatial quantile regression
#'
#' Local bivariate spline estimation, decorrelated initialization, and
#' communication for quantile spatial models with constant and varying terms.
#'
#' @details
#' Use \code{\link{subregion_const}} to form worker subregions and
#' \code{\link{forward.disq}} to compare coefficient allocations.
#' \code{\link{disq_from_local_fits}} combines saved local SQBiT estimates
#' with a fixed global projection and provides spatial covariance intervals.
#' \code{\link{disq_main}} retains the earlier fitting interface.
#'
#' A simulated example is installed at
#' \code{system.file("examples", "getting-started.R", package = "DISQ")}.
#'
#' @name DISQ-package
#' @aliases DISQ
"_PACKAGE"
