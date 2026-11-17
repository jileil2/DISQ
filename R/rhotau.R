#' Quantile check loss
#'
#' Evaluates \eqn{u\{\tau-I(u<0)\}} for each residual.
#' @param u Numeric residual vector.
#' @param tau Quantile level.
#' @return The check loss for each element of \code{u}.
#' @export
rhotau <- function(u, tau) {
  return(u * (tau - (u < 0)))
}
