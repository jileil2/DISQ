# Kernel scores, local weighted projection, and root solvers.
.disq_kernel_inputs <- function(value, h, kernel) {
  kernel <- match.arg(kernel, c("unif", "norm"))
  if (!is.numeric(value) || any(!is.finite(value)) ||
    length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("Kernel inputs must be finite with a positive bandwidth.")
  }
  kernel
}

.disq_kernel_density <- function(residual, h, kernel = "unif") {
  kernel <- .disq_kernel_inputs(residual, h, kernel)
  if (kernel == "unif") dunif(residual / h, -1, 1) / h else
    dnorm(residual / h) / h
}

.disq_kernel_cdf <- function(value, h, kernel = "unif") {
  kernel <- .disq_kernel_inputs(value, h, kernel)
  if (kernel == "unif") punif(value / h, -1, 1) else pnorm(value / h)
}

.disq_local_projection <- function(Z, C, residual, h, weights = NULL,
                                   kernel = "unif") {
  Z <- as.matrix(Z)
  C <- as.matrix(C)
  residual <- as.numeric(residual)
  n <- nrow(Z)
  L <- ncol(Z)
  if (n < 1L || L < 1L || ncol(C) < 1L || nrow(C) != n ||
    length(residual) != n || !is.numeric(Z) || !is.numeric(C) ||
    any(!is.finite(Z)) || any(!is.finite(C)) ||
    any(!is.finite(residual))) {
    stop("Z, C, and residual must be finite and have matching observations.")
  }
  if (length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("h must be a finite positive bandwidth.")
  }
  kernel <- match.arg(kernel, c("unif", "norm"))
  if (is.null(weights)) {
    weights <- .disq_kernel_density(residual, h, kernel)
  }
  weights <- as.numeric(weights)
  if (length(weights) != n || any(!is.finite(weights)) ||
    any(weights < 0)) {
    stop("Projection weights must be finite, nonnegative, and length n.")
  }
  positive <- which(weights > 0)
  rank_error <- function(rank) {
    stop(sprintf(
      paste0(
        "The weighted spline design is singular or numerically rank deficient: ",
        "rank %d/%d with %d positive weights; ",
        "the full-column-rank projection is undefined."
      ),
      rank, L, length(positive)
    ))
  }
  if (!length(positive)) rank_error(0L)
  weighted_design <- Z[positive, , drop = FALSE] * sqrt(weights[positive])
  weighted_covariates <- C[positive, , drop = FALSE] * sqrt(weights[positive])
  scales <- sqrt(colSums(weighted_design^2))
  if (any(!is.finite(scales))) {
    stop("Weighted spline column norms are not finite.")
  }
  # Zero columns stay zero during rank detection. LAPACK's qr$rank is not a
  # numerical-rank diagnostic, so inspect the pivoted R diagonals explicitly.
  scales[scales == 0] <- 1
  decomposition <- qr(sweep(weighted_design, 2L, scales, "/"), LAPACK = TRUE)
  factor <- qr.R(decomposition)
  diagonals <- abs(diag(factor))
  rank_tolerance <- max(dim(weighted_design)) * .Machine$double.eps *
    max(diagonals)
  rank <- sum(diagonals > rank_tolerance)
  if (rank != L) rank_error(rank)
  coefficients <- sweep(
    qr.coef(decomposition, weighted_covariates),
    1L, scales, "/"
  )
  residualized <- C - Z %*% coefficients
  rhs <- crossprod(Z, C * weights) / n
  normal_error <- crossprod(Z, residualized * weights) / n
  orthogonality_error <- max(abs(normal_error)) / max(1, max(abs(rhs)))

  list(
    P = t(coefficients),
    residualized = residualized,
    weights = weights,
    orthogonality_error = orthogonality_error,
    rank = rank,
    # Undo column scaling in R. Its squared condition estimates the condition
    # of the original weighted Gram matrix, as in the earlier implementation.
    condition = kappa(sweep(factor, 2L, scales[decomposition$pivot], "*"),
      exact = FALSE, method = "direct", triangular = TRUE
    )^2
  )
}

.disq_score_inputs <- function(eta, Y, C, nuisance, residualized, h, tau) {
  C <- as.matrix(C)
  residualized <- as.matrix(residualized)
  eta <- as.numeric(eta)
  Y <- as.numeric(Y)
  nuisance <- as.numeric(nuisance)
  n <- nrow(C)
  q <- ncol(C)
  if (n < 1L || q < 1L || length(Y) != n || length(nuisance) != n ||
    length(eta) != q || !identical(dim(residualized), dim(C)) ||
    !is.numeric(C) || !is.numeric(residualized) ||
    any(!is.finite(C)) || any(!is.finite(residualized)) ||
    any(!is.finite(Y)) || any(!is.finite(nuisance)) ||
    any(!is.finite(eta))) {
    stop("Score inputs must be finite with matching n observations and q coefficients.")
  }
  if (length(h) != 1L || !is.finite(h) || h <= 0 ||
    length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("Use a positive bandwidth and a quantile level strictly between 0 and 1.")
  }
  list(
    eta = eta, Y = Y, C = C, nuisance = nuisance,
    residualized = residualized, h = h, tau = tau, n = n, q = q
  )
}

.disq_score <- function(eta, Y, C, nuisance, residualized, h, tau,
                        kernel = "unif") {
  dat <- .disq_score_inputs(eta, Y, C, nuisance, residualized, h, tau)
  residual <- dat$Y - dat$nuisance - as.numeric(dat$C %*% dat$eta)
  psi <- .disq_kernel_cdf(-residual, h, kernel) - tau
  as.numeric(colMeans(dat$residualized * psi))
}

.disq_jacobian <- function(eta, Y, C, nuisance, residualized, h, tau,
                           kernel = "unif") {
  dat <- .disq_score_inputs(eta, Y, C, nuisance, residualized, h, tau)
  residual <- dat$Y - dat$nuisance - as.numeric(dat$C %*% dat$eta)
  weights <- .disq_kernel_density(residual, h, kernel)
  crossprod(dat$residualized, dat$C * weights) / dat$n
}

# Solve the fixed-instrument score equation.  The master Jacobian need not be
# symmetric, so the Newton direction uses an ordinary linear-system solve.
.disq_solve_score <- function(initial, Y, C, nuisance, residualized,
                              h, tau, offset = NULL, max.iter = 100,
                              tolerance = 1e-8, kernel = "unif") {
  dat <- .disq_score_inputs(initial, Y, C, nuisance, residualized, h, tau)
  kernel <- match.arg(kernel, c("unif", "norm"))
  if (length(max.iter) != 1L || !is.finite(max.iter) || max.iter < 1 ||
    max.iter != as.integer(max.iter) || length(tolerance) != 1L ||
    !is.finite(tolerance) || tolerance <= 0) {
    stop("max.iter must be a positive integer and tolerance must be positive.")
  }
  if (is.null(offset)) offset <- numeric(dat$q)
  offset <- as.numeric(offset)
  if (length(offset) != dat$q || any(!is.finite(offset))) {
    stop("offset must be a finite vector of length q.")
  }
  score <- function(eta) {
    residual <- dat$Y - dat$nuisance - as.numeric(dat$C %*% eta)
    psi <- .disq_kernel_cdf(-residual, h, kernel) - tau
    as.numeric(colMeans(dat$residualized * psi)) - offset
  }
  jacobian <- function(eta) {
    residual <- dat$Y - dat$nuisance - as.numeric(dat$C %*% eta)
    weights <- .disq_kernel_density(residual, h, kernel)
    crossprod(dat$residualized, dat$C * weights) / dat$n
  }

  eta <- dat$eta
  value <- score(eta)
  iterations <- 0L
  method <- "damped Newton"
  status <- "maximum iterations reached"
  stalled <- FALSE
  while (max(abs(value)) > tolerance && iterations < max.iter) {
    direction <- tryCatch(as.numeric(solve(jacobian(eta), value)),
      error = function(e) NULL
    )
    if (is.null(direction) || any(!is.finite(direction))) {
      status <- "singular Jacobian"
      stalled <- TRUE
      break
    }
    merit <- sum(value^2)
    accepted <- FALSE
    damping <- 1
    for (backtrack in seq_len(30L)) {
      candidate <- eta - damping * direction
      if (all(is.finite(candidate))) {
        candidate_value <- score(candidate)
        if (all(is.finite(candidate_value)) &&
          (max(abs(candidate_value)) <= tolerance ||
            sum(candidate_value^2) <= (1 - 1e-4 * damping) * merit)) {
          accepted <- TRUE
          break
        }
      }
      damping <- damping / 2
    }
    if (!accepted) {
      status <- "Newton line search stalled"
      stalled <- TRUE
      break
    }
    eta <- candidate
    value <- candidate_value
    iterations <- iterations + 1L
  }

  if (stalled && max(abs(value)) > tolerance && iterations < max.iter &&
    requireNamespace("nleqslv", quietly = TRUE)) {
    fallback <- tryCatch(
      nleqslv::nleqslv(
        eta, score,
        jac = jacobian, method = "Newton", global = "dbldog",
        control = list(
          ftol = tolerance, xtol = .Machine$double.eps^(2 / 3),
          maxit = as.integer(max.iter - iterations)
        )
      ),
      error = function(e) NULL
    )
    if (!is.null(fallback) && all(is.finite(fallback$x))) {
      fallback_value <- score(fallback$x)
      iterations <- iterations + fallback$iter
      if (max(abs(fallback_value)) < max(abs(value))) {
        eta <- as.numeric(fallback$x)
        value <- fallback_value
      }
      method <- "damped Newton / nleqslv"
      status <- fallback$message
    }
  }

  score_norm <- max(abs(value))
  converged <- is.finite(score_norm) && score_norm <= tolerance
  if (converged) status <- "score tolerance reached"
  list(
    eta = eta, converged = converged, score_norm = score_norm,
    iterations = iterations, score = value, method = method, status = status
  )
}
