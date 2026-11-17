#' Generate quantile wild-bootstrap weights
#'
#' Draws weights equal to \eqn{-2\tau} with probability \eqn{\tau} and
#' \eqn{2(1-\tau)} otherwise.
#' @param n Number of weights.
#' @param tau Quantile level.
#' @return A numeric vector of length \code{n}.
#' @export
rG <- function(n, tau) {
  Is <- rbinom(n, size = 1, prob = tau)
  Is * (-2 * tau) + 2 * (1 - tau) * (1 - Is)
}
