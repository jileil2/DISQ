# ============================================================
# smqr.eta.dc(): Barzilai-Borwein optimizer for decorrelated
#                smoothed quantile regression
# ============================================================

#' Barzilai-Borwein Optimization for Decorrelated Smoothed Quantile Regression
#'
#' \code{smqr.eta.dc()} updates the constant-coefficient vector in a
#' decorrelated smoothed quantile regression problem using a gradient-based
#' algorithm with Barzilai-Borwein step-size selection.
#'
#' This routine is used in the communication step of DISQ to minimize the
#' surrogate objective for the constant coefficients. The objective function
#' is not formed explicitly. Instead, the algorithm iteratively updates the
#' parameter vector using the gradient of the decorrelated smoothed score.
#'
#' Let \eqn{\eta} denote the constant-coefficient vector. For a given
#' bandwidth \eqn{h}, define the pseudo-response
#' \deqn{
#' Y^* = Y - W\gamma,
#' }
#' the residual vector
#' \deqn{
#' r = Y^* - C\eta,
#' }
#' and the projected constant-covariate matrix
#' \deqn{
#' C^* = GC.
#' }
#' The gradient used by the algorithm is
#' \deqn{
#' \nabla(\eta)
#' =
#' \frac{1}{n} (C^*)^\top
#' \left\{
#' F_U(-r/h) - \tau
#' \right\}
#' - \mathrm{grad},
#' }
#' where \eqn{F_U(\cdot)} is the cumulative distribution function of the
#' standard uniform distribution, corresponding to the integrated uniform
#' kernel used in the convolution-smoothed quantile loss.
#'
#' Starting from \code{initial}, the algorithm performs gradient updates
#' with a Barzilai-Borwein step size computed from successive differences in
#' the parameter and gradient vectors. This provides a simple curvature-based
#' acceleration without explicitly forming a Hessian matrix.
#'
#' @param h Positive smoothing bandwidth used in the convolution-smoothed
#'   quantile loss.
#' @param tau Quantile level in \eqn{(0,1)}.
#' @param Y Numeric response vector of length \eqn{n}.
#' @param C Numeric design matrix for the constant coefficients.
#' @param W Numeric design matrix for the varying-coefficient component.
#' @param G Numeric projection matrix used in the decorrelated score
#'   construction.
#' @param grad Numeric gradient offset vector. In DISQ, this is typically
#'   the difference between the master-machine gradient and the aggregated
#'   global gradient.
#' @param initial Initial value of the constant-coefficient vector.
#' @param gamma Numeric vector of nuisance parameters associated with the
#'   varying-coefficient component.
#' @param max.iter Maximum number of gradient iterations.
#' @param stepMax Upper bound for the Barzilai-Borwein step size.
#' @param epsilon Convergence tolerance based on the Frobenius norm of the
#'   parameter update.
#'
#' @return A column matrix containing the updated estimate of the constant
#'   coefficients.
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item compute the decorrelated smoothed gradient at the current iterate;
#'   \item take one gradient step;
#'   \item update the step size using the Barzilai-Borwein rule;
#'   \item repeat until the update norm is below \code{epsilon} or the
#'   iteration limit is reached.
#' }
#'
#' The Barzilai-Borwein step size is computed from
#' \eqn{\eta^{(t)} - \eta^{(t-1)}} and
#' \eqn{\nabla^{(t)} - \nabla^{(t-1)}}, which approximates second-order
#' curvature information while remaining computationally light.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' n <- 200
#' q <- 3
#' p <- 5
#'
#' Y <- rnorm(n)
#' C <- matrix(rnorm(n * q), nrow = n, ncol = q)
#' W <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' G <- diag(n)
#' gamma <- matrix(0, nrow = p, ncol = 1)
#' grad <- matrix(0, nrow = q, ncol = 1)
#' initial <- matrix(0, nrow = q, ncol = 1)
#'
#' eta_hat <- smqr.eta.dc(
#'   h = 0.5,
#'   tau = 0.5,
#'   Y = Y,
#'   C = C,
#'   W = W,
#'   G = G,
#'   grad = grad,
#'   initial = initial,
#'   gamma = gamma,
#'   max.iter = 100,
#'   stepMax = 20,
#'   epsilon = 1e-5
#' )
#' }
#'
#' @export
smqr.eta.dc <- function(h, tau, Y, C, W, G, grad, initial, gamma,
                        max.iter = 100, stepMax = 20, epsilon = 1e-5) {
  # ---------------------------
  # Basic checks and coercions
  # ---------------------------
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("smqr.eta.dc(): 'h' must be a positive finite scalar.")
  }
  if (!is.numeric(tau) || length(tau) != 1L || tau <= 0 || tau >= 1) {
    stop("smqr.eta.dc(): 'tau' must be a scalar in (0, 1).")
  }
  if (!is.numeric(max.iter) || length(max.iter) != 1L || max.iter < 1) {
    stop("smqr.eta.dc(): 'max.iter' must be a positive integer.")
  }
  if (!is.numeric(stepMax) || length(stepMax) != 1L || !is.finite(stepMax) || stepMax <= 0) {
    stop("smqr.eta.dc(): 'stepMax' must be a positive finite scalar.")
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
    stop("smqr.eta.dc(): 'epsilon' must be a positive finite scalar.")
  }

  Y <- as.numeric(Y)
  C <- as.matrix(C)
  W <- as.matrix(W)
  G <- as.matrix(G)
  grad <- matrix(as.numeric(grad), ncol = 1)
  initial <- matrix(as.numeric(initial), ncol = 1)
  gamma <- matrix(as.numeric(gamma), ncol = 1)

  n <- length(Y)

  if (nrow(C) != n) {
    stop("smqr.eta.dc(): nrow(C) must equal length(Y).")
  }
  if (nrow(W) != n) {
    stop("smqr.eta.dc(): nrow(W) must equal length(Y).")
  }
  if (nrow(G) != n || ncol(G) != n) {
    stop("smqr.eta.dc(): G must be an n x n matrix, where n = length(Y).")
  }
  if (ncol(C) != nrow(initial)) {
    stop("smqr.eta.dc(): length of 'initial' must match ncol(C).")
  }
  if (ncol(C) != nrow(grad)) {
    stop("smqr.eta.dc(): length of 'grad' must match ncol(C).")
  }
  if (ncol(W) != nrow(gamma)) {
    stop("smqr.eta.dc(): length of 'gamma' must match ncol(W).")
  }

  # ---------------------------
  # Gradient function
  # ---------------------------
  nabla <- function(eta, h) {
    Y.star <- Y - as.numeric(eigenMapMatMult(W, gamma))
    r <- Y.star - as.numeric(eigenMapMatMult(C, eta))
    C.star <- as.matrix(eigenMapMatMult(G, C))
    sc <- matrix(punif(-r / h, -1, 1) - tau, ncol = 1)
    as.matrix(eigenMapMatMult(t(C.star), sc)) / n - grad
  }

  # ---------------------------
  # Initialization
  # ---------------------------
  eta <- initial
  grad_old <- nabla(eta = eta, h = h)

  eta <- eta - grad_old
  grad_new <- nabla(eta = eta, h = h)

  eta_diff <- -grad_old
  grad_diff <- grad_new - grad_old
  iter <- 1L

  # ---------------------------
  # Barzilai-Borwein iterations
  # ---------------------------
  while (base::norm(eta_diff, type = "F") > epsilon && iter <= max.iter) {
    alpha <- 1

    cross <- as.numeric(t(eta_diff) %*% grad_diff)
    grad_norm_sq <- as.numeric(t(grad_diff) %*% grad_diff)
    eta_norm_sq <- as.numeric(t(eta_diff) %*% eta_diff)

    if (is.finite(cross) && is.finite(grad_norm_sq) &&
        cross > 0 && grad_norm_sq > 0) {
      a1 <- cross / grad_norm_sq
      a2 <- eta_norm_sq / cross
      alpha <- min(c(a1, a2, stepMax))
    }

    grad_old <- grad_new
    eta_diff <- -alpha * grad_new
    eta <- eta + eta_diff
    grad_new <- nabla(eta = eta, h = h)
    grad_diff <- grad_new - grad_old
    iter <- iter + 1L
  }

  eta
}
