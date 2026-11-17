# Finite-range spatial covariance for the constant coefficients.
.disq_pair_kernel <- function(S, radius) {
  S <- as.matrix(S)
  n <- nrow(S)
  if (n < 1L || ncol(S) != 2L || any(!is.finite(S)) ||
    length(radius) != 1L || !is.finite(radius) || radius <= 0) {
    stop("A finite two-dimensional spatial design and positive cutoff are required.")
  }
  if (!requireNamespace("dbscan", quietly = TRUE) ||
    !requireNamespace("Matrix", quietly = TRUE)) {
    stop("The finite-range covariance kernel requires dbscan and Matrix.")
  }
  neighbors <- dbscan::frNN(S, eps = radius, sort = FALSE)$id
  neighbors <- Map(function(ids, i) unique(c(i, ids)), neighbors, seq_len(n))
  Matrix::sparseMatrix(
    i = rep(seq_len(n), lengths(neighbors)),
    j = as.integer(unlist(neighbors, use.names = FALSE)),
    x = 1, dims = c(n, n)
  )
}

.disq_covariance <- function(eta, Y, C, nuisance, D, h, tau, pair_kernel,
                             level = 0.05, kernel = "norm") {
  n <- length(Y)
  q <- ncol(C)
  residual <- as.numeric(Y - nuisance - C %*% eta)
  score <- .disq_kernel_cdf(-residual, h, kernel) - tau
  U <- D * score
  J <- crossprod(D, C * .disq_kernel_density(residual, h, kernel)) / n
  Gamma <- as.matrix(crossprod(U, pair_kernel %*% U)) / n
  unavailable <- matrix(NA_real_, q, q)
  V <- unavailable
  se <- rep(NA_real_, q)
  ci <- matrix(NA_real_, q, 2L)
  inverse <- tryCatch(solve(J), error = function(e) e)
  problem <- NULL
  if (inherits(inverse, "error")) {
    problem <- paste("Global Jacobian solve failed:", conditionMessage(inverse))
  } else {
    V <- inverse %*% Gamma %*% t(inverse) / n
    # Remove only floating-point antisymmetry, never truncate eigenvalues.
    V <- (V + t(V)) / 2
    if (any(!is.finite(V))) {
      problem <- "Estimated covariance contains nonfinite entries."
    } else {
      minimum_eigenvalue <- min(eigen(V, symmetric = TRUE, only.values = TRUE)$values)
      if (any(diag(V) <= 0) || minimum_eigenvalue <= 0) {
        problem <- "Estimated covariance is not positive definite."
      } else {
        se <- sqrt(diag(V))
        critical <- stats::qnorm(1 - level / 2)
        ci <- cbind(eta - critical * se, eta + critical * se)
      }
    }
  }
  list(
    V = V, J = J, Gamma = Gamma, se = se, ci = ci,
    valid = is.null(problem), error = problem,
    diagnostics = list(
      jacobian_rcond = tryCatch(rcond(J), error = function(e) NA_real_),
      covariance_min_eigenvalue = if (all(is.finite(V)))
        min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) else NA_real_,
      gamma_min_eigenvalue = if (all(is.finite(Gamma)))
        min(eigen((Gamma + t(Gamma)) / 2,
          symmetric = TRUE,
          only.values = TRUE
        )$values) else NA_real_
    )
  )
}
