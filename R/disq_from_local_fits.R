disq_from_local_fits <- function(Y, C, X, S, B, Q2, index, local, tau,
                                 radius, master = which.max(lengths(index)),
                                 c = 0.20, rounds = 15L, level = 0.05,
                                 score_tolerance = 1e-7, B_eval = NULL) {
  Y <- as.numeric(Y)
  C <- as.matrix(C)
  X <- as.matrix(X)
  S <- as.matrix(S)
  Q2 <- as.matrix(Q2)
  N <- length(Y)
  ns <- lengths(index)
  ids <- unlist(index, use.names = FALSE)
  if (N < 1L || !is.numeric(ids) ||
    !identical(as.numeric(sort(ids, na.last = TRUE)), as.numeric(seq_len(N))) ||
    length(local) != length(index) || any(ns == 0L)) {
    stop("index must partition the observations into nonempty workers.")
  }
  if (nrow(C) != N || nrow(X) != N || nrow(S) != N || ncol(S) != 2L ||
    nrow(B) != N || nrow(Q2) != ncol(B) || ncol(Q2) < 1L || ncol(X) < 1L ||
    any(!is.finite(Y)) || any(!is.finite(C)) || any(!is.finite(X)) ||
    any(!is.finite(S)) || !is.finite(max(abs(B))) || any(!is.finite(Q2))) {
    stop("Supply finite inputs with matching observation and spline dimensions.")
  }
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("tau must lie strictly between zero and one.")
  }

  raw <- matrix(0, ncol(B), ncol(X))
  columns <- lapply(local, function(fit) fit$columns)
  column_ids <- unlist(columns, use.names = FALSE)
  if (!is.numeric(column_ids) ||
    !identical(as.numeric(sort(column_ids, na.last = TRUE)), as.numeric(seq_len(ncol(B))))) {
    stop("Local raw spline columns must partition the global raw spline basis.")
  }
  for (k in seq_along(local)) {
    theta <- as.matrix(local[[k]]$theta)
    if (!identical(dim(theta), c(length(columns[[k]]), ncol(X))) ||
      any(!is.finite(theta))) {
      stop("Invalid raw local spline coefficients on worker ", k, ".")
    }
    raw[columns[[k]], ] <- theta
  }
  corrected <- Q2 %*% crossprod(Q2, raw)
  varying <- as.matrix(B %*% corrected)
  nuisance <- rowSums(X * varying)
  evaluated <- if (is.null(B_eval)) NULL else as.matrix(B_eval %*% corrected)
  if (ncol(C) == 0L) {
    return(list(
      eta = numeric(0), theta = corrected, beta = varying,
      beta_eval = evaluated, fitted = nuisance, completed = TRUE
    ))
  }
  if (length(master) != 1L || !is.finite(master) || master != floor(master) ||
    master < 1L || master > length(index)) {
    stop("Choose a valid master worker.")
  }

  local_bases <- lapply(seq_along(local), function(k) {
    basis <- B[index[[k]], columns[[k]], drop = FALSE]
    constraint <- as.matrix(local[[k]]$Q2)
    if (nrow(constraint) != ncol(basis) || ncol(constraint) < 1L ||
      any(!is.finite(constraint))) {
      stop("Invalid local constrained basis on worker ", k, ".")
    }
    as.matrix(basis %*% constraint)
  })
  local_designs <- lapply(seq_along(local), function(k) {
    x <- X[index[[k]], , drop = FALSE]
    do.call(cbind, lapply(seq_len(ncol(x)), function(j) local_bases[[k]] * x[, j]))
  })
  local_nuisance <- lapply(seq_along(local), function(k) {
    beta <- B[index[[k]], columns[[k]], drop = FALSE] %*% raw[columns[[k]], , drop = FALSE]
    rowSums(X[index[[k]], , drop = FALSE] * as.matrix(beta))
  })
  eta_local <- lapply(local, function(fit) as.numeric(fit$eta))
  if (any(lengths(eta_local) != ncol(C))) {
    stop("Each local fit must have one estimate per column of C.")
  }
  eta_local <- do.call(rbind, eta_local)
  phi_local <- vapply(local, function(fit) fit$phi, numeric(1))
  initialization <- .disq_local_roots(
    lapply(index, function(i) Y[i]), lapply(index, function(i) C[i, , drop = FALSE]),
    local_designs, local_nuisance, eta_local, phi_local, tau
  )
  bw <- disq_bandwidths(N, ncol(Q2), ns[master], ncol(local_bases[[master]]), c)
  weighted_residual <- Y - nuisance - as.numeric(C %*% initialization$eta_mean)
  projection <- .disq_global_projection(
    B, Q2, X, C, .disq_kernel_density(weighted_residual, bw[["h"]], "norm")
  )
  iterates <- .disq_communicate(
    Y, C, S, nuisance, projection$D, index, ns, master,
    bw[["h"]], bw[["b"]], tau,
    max_rounds = rounds, tol = score_tolerance,
    vcov_radius = radius, eta_initial = initialization$eta_mean, level = level
  )
  final_eta <- if (iterates$completed) iterates$eta else rep(NA_real_, ncol(C))
  list(
    eta = final_eta, theta = corrected, beta = varying, beta_eval = evaluated,
    fitted = nuisance + as.numeric(C %*% final_eta),
    completed = iterates$completed, converged = iterates$converged,
    bandwidths = bw, initialization = initialization, projection = projection,
    communication = iterates,
    covariance = if (iterates$completed) tail(iterates$covariance, 1L)[[1L]] else NULL
  )
}
