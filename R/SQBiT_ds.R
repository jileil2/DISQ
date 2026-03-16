# ============================================================
# SQBiT_ds(): patched SQBiT_ds to avoid "Not a matrix."
# ============================================================

#' Decorrelated-Score ADMM Solver (Safe / Type-Stable Version)
#'
#' \code{SQBiT_ds()} is a **patched, type-stable** variant of the internal DISQ routine
#' used to solve the decorrelated-score estimating equation for the constant effects
#' after subtracting the estimated varying-coefficient component.
#'
#' The original implementation can fail on some platforms / parallel backends with
#' errors such as \dQuote{Not a matrix.} due to implicit type conversions
#' (e.g., numeric vectors becoming \code{NULL}, 1D objects, or integer storage modes).
#' This safe variant enforces:
#' \itemize{
#'   \item coercion of \code{C}, \code{X}, \code{B}, \code{Q2} to \code{double} matrices,
#'   \item explicit matrix shapes for intermediate objects (e.g., \code{rc}, \code{gamma.0}),
#'   \item robust matrix multiplication via \code{eigenMapMatMult()} with fallback to \code{\%*\%}.
#' }
#'
#' Internally, it forms
#' \eqn{Z = \mathrm{kr}(X, BQ_2)} and \eqn{W = [C, Z]},
#' then computes the weighted projection matrix \eqn{G} (onto the orthogonal complement
#' of the varying-coefficient design, with kernel weights),
#' and finally fits a weighted, convolution-smoothed quantile regression via ADMM.
#'
#' @param tau Quantile level in \eqn{(0,1)}. If outside \eqn{[0,1]}, the function resets \code{tau = 0.5}
#' and issues a warning.
#' @param Y Numeric response vector of length \eqn{n}.
#' @param C Numeric matrix (\eqn{n \times q}) of covariates with constant effects.
#' @param X Numeric matrix (\eqn{n \times p}) of covariates with varying coefficients.
#' @param B Bernstein basis matrix (\eqn{n \times J_B}) evaluated at observation locations.
#' @param Q2 Matrix \eqn{Q_2} from the QR decomposition used to enforce spline smoothness.
#' @param h Bandwidth for convolution smoothing.
#' @param lambda Penalization parameter(s). Accepts a scalar, vector, or matrix; it is expanded
#' internally to an \eqn{n_\lambda \times p} matrix if needed.
#' @param eta.0 Initial value for constant coefficients \eqn{\eta} (length \eqn{q}).
#' @param gamma.0 Initial value for spline coefficients \eqn{\gamma} (vector).
#' @param s.j ADMM step size parameter (a.k.a. \eqn{\rho}). Default is \code{1}.
#' @param eps.abs Absolute tolerance for ADMM stopping criterion.
#' @param eps.rel Relative tolerance for ADMM stopping criterion.
#' @param var.j Logical; if \code{TRUE}, adaptively updates \code{s.j} using \code{zeta} and \code{incr}.
#' @param zeta Threshold parameter for adaptive step size updates.
#' @param incr Multiplicative factor for step size updates.
#' @param max.iter Maximum ADMM iterations.
#'
#' @details
#' This function assumes the following compiled / package-level helpers exist in the calling environment:
#' \code{eigenMapMatMult()} and \code{armaInv()}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{eta}}{Estimated constant-effect vector (including intercept if present in \code{C}).}
#'   \item{\code{iter}}{Number of ADMM iterations used.}
#'   \item{\code{G}}{Weighted projection matrix used by the decorrelation step.}
#'   \item{\code{time}}{Elapsed runtime.}
#' }
#'
#' @keywords internal
SQBiT_ds <- function(tau = .5, Y, C, X, B, Q2, h, lambda,
                          eta.0, gamma.0, s.j = 1,
                          eps.abs = 1e-4, eps.rel = 1e-2,
                          var.j = FALSE, zeta = 10, incr = 2,
                          max.iter = 50) {

  start_time <- Sys.time()

  Y  <- as.numeric(Y)
  C  <- as.matrix(C);  storage.mode(C)  <- "double"
  X  <- as.matrix(X);  storage.mode(X)  <- "double"
  B  <- as.matrix(B);  storage.mode(B)  <- "double"
  Q2 <- as.matrix(Q2); storage.mode(Q2) <- "double"

  eta.0   <- as.numeric(eta.0)
  gamma.0 <- as.numeric(gamma.0)

  eMM <- function(A, B) {
    A <- as.matrix(A); storage.mode(A) <- "double"
    B <- as.matrix(B); storage.mode(B) <- "double"
    out <- try(eigenMapMatMult(A, B), silent = TRUE)
    if (inherits(out, "try-error")) out <- A %*% B
    out
  }

  if (tau < 0 || tau > 1) {
    warning("Quantile level should be between 0 and 1. Using tau=0.5.")
    tau <- 0.5
  }

  n  <- length(Y)
  np <- ncol(X)
  nq <- ncol(C)

  lambda <- as.matrix(lambda); storage.mode(lambda) <- "double"
  nl <- nrow(lambda)
  if (ncol(lambda) == 1) lambda <- matrix(rep(lambda, times = np), nl, np)

  BQ2 <- eMM(B, Q2)
  Z <- as.matrix(kr(X, BQ2, byrow = TRUE)); storage.mode(Z) <- "double"
  W <- as.matrix(cbind(C, Z));            storage.mode(W) <- "double"

  rc <- matrix(c(eta.0, gamma.0), ncol = 1); storage.mode(rc) <- "double"

  rs <- Y - as.numeric(eMM(W, rc))
  kd <- as.numeric(dunif(rs / h, -1, 1) / h)

  W.h.st <- as.matrix(diag(sqrt(kd))); storage.mode(W.h.st) <- "double"
  st.X <- eMM(W.h.st, X)

  X.Wh.X <- eMM(t(st.X), st.X)
  X.Wh.X.inv <- armaInv(X.Wh.X)

  G <- eMM(X.Wh.X.inv, t(st.X))
  G <- diag(n) - eMM(st.X, G)
  G <- as.matrix(G); storage.mode(G) <- "double"

  ws <- as.numeric(rowSums(G))
  ws <- ws / n

  Y.star <- Y - as.numeric(eMM(Z, matrix(gamma.0, ncol = 1)))

  C1 <- as.matrix(cbind(1, C)); storage.mode(C1) <- "double"
  mat <- armaInv(eMM(t(C1), C1))
  mat <- eMM(mat, t(C1))

  psi.eq <- function(r, s.j) {
    psi1 <- (s.j * r - (tau - 1 / 2) * ws) / (s.j + ws / (2 * h))
    psi2 <- - (tau - 1) / s.j * ws + r
    psi3 <- - tau / s.j * ws + r
    psi1.r <- psi1 * (abs(psi1) <= h)
    psi2.r <- psi2 * (psi2 < -h)
    psi3.r <- psi3 * (psi3 >  h)
    psi1.r + psi2.r + psi3.r
  }

  new.u <- (1 / s.j) * Y.star
  new.eta <- matrix(c(0, eta.0), ncol = 1)
  new.fit <- as.numeric(eMM(C1, new.eta))

  iter.admm <- 0
  repeat {
    iter.admm <- iter.admm + 1

    old.u <- new.u
    old.fit <- new.fit

    r <- Y.star - old.fit - old.u
    new.psi <- psi.eq(r = r, s.j = s.j)

    new.eta <- eMM(mat, matrix(Y.star - new.psi - old.u, ncol = 1))
    new.fit <- as.numeric(eMM(C1, new.eta))
    res.new <- Y.star - new.fit
    new.u <- old.u + (new.psi - res.new)

    r.prime <- res.new - new.psi
    r.dual <- s.j * (new.fit - old.fit)

    eps.prime <- sqrt(length(r.prime)) * eps.abs +
      eps.rel * max(norm(new.psi, type = "2"),
                    norm(new.fit, type = "2"),
                    norm(Y.star, type = "2"))
    eps.dual <- sqrt(length(r.dual)) * eps.abs +
      eps.rel * norm(new.u, type = "2")

    n.prime <- norm(r.prime, type = "2")
    n.dual  <- norm(r.dual,  type = "2")

    if (isTRUE(var.j)) {
      old.s.j <- s.j
      if (n.prime > zeta * n.dual) s.j <- incr * old.s.j
      else if (zeta * n.prime < n.dual) s.j <- old.s.j / incr
      else s.j <- old.s.j
    }

    if ((n.prime < eps.prime) && (n.dual < eps.dual)) break
    if (iter.admm >= max.iter) break
  }

  end_time <- Sys.time()

  list(eta = new.eta, iter = iter.admm, G = G, time = end_time - start_time)
}
