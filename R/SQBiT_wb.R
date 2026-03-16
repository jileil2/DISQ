# ============================================================
# SQBiT_wb(): smoothed QSVCM with wild bootstrap
# ============================================================

#' Wild-bootstrap inference for smoothed quantile spatial varying coefficient models
#'
#' Fits a smoothed quantile spatial varying coefficient model (QSVCM) and
#' performs wild-bootstrap inference for both constant coefficients and
#' varying-coefficient functions.
#'
#' @param h A positive numeric bandwidth used in the smoothed quantile loss.
#' @param tau A numeric scalar in (0,1) specifying the quantile level.
#' @param Y A numeric response vector of length n.
#' @param C A design matrix for constant-effect covariates with n rows.
#' @param X A design matrix for varying-effect covariates with n rows.
#' @param P A penalty matrix for the spline coefficients.
#' @param B A spline basis matrix evaluated at the observed spatial locations.
#' @param Q2 A reduced-basis or constraint matrix used in the spline representation.
#' @param max.iter Maximum number of ADMM iterations for each bootstrap fit.
#' @param lambda A smoothing parameter object. It is sanitized internally to a
#'   numeric vector of length ncol(X).
#' @param eps.abs Absolute tolerance for the ADMM stopping criterion.
#' @param eps.rel Relative tolerance for the ADMM stopping criterion.
#' @param kernel Character string specifying the smoothing kernel. Supported
#'   values are "unif" and "norm".
#' @param var.j Logical; if TRUE, adaptively update eta.j.
#' @param zeta Positive numeric tuning parameter controlling adaptive eta.j updates.
#' @param incr Multiplicative adjustment factor for adaptive eta.j updates.
#' @param eta.j Initial ADMM penalty parameter.
#' @param BQ2.eva Basis evaluation matrix for varying coefficients.
#' @param compute.vc Logical; if TRUE, compute bootstrap inference for varying coefficients.
#' @param Br Number of wild-bootstrap replications.
#' @param level Significance level for confidence intervals.
#' @param gamma.hat Fitted coefficient vector from the original sample.
#' @param biascorr Logical; if TRUE, apply bias correction.
#'
#' @return A list containing bootstrap draws and confidence intervals.
#' @export
SQBiT_wb <- function(h, tau, Y, C, X, P, B, Q2, max.iter = 50,
                     lambda, eps.abs = 1e-4, eps.rel = 1e-2,
                     kernel = "unif", var.j = FALSE, zeta = 10,
                     incr = 2, eta.j = 0.32, BQ2.eva,
                     compute.vc = TRUE, Br = 200,
                     level = 0.05, gamma.hat, biascorr = FALSE) {
  start_time <- Sys.time()

  Y  <- as.numeric(Y)
  C  <- as.matrix(C)
  X  <- as.matrix(X)
  P  <- as.matrix(P)
  B  <- as.matrix(B)
  Q2 <- as.matrix(Q2)
  BQ2.eva <- as.matrix(BQ2.eva)

  if (tau <= 0 || tau >= 1) {
    warning("Quantile level should be between 0 and 1. The default quantile level of 0.5 is used.")
    tau <- 0.5
  }

  n  <- length(Y)
  np <- ncol(X)
  nq <- ncol(C)
  J  <- ncol(Q2)

  if (nrow(C) != n) stop("SQBiT_wb(): nrow(C) must equal length(Y).")
  if (nrow(X) != n) stop("SQBiT_wb(): nrow(X) must equal length(Y).")
  if (nrow(B) != n) stop("SQBiT_wb(): nrow(B) must equal length(Y).")
  if (nrow(P) != J || ncol(P) != J) {
    stop("SQBiT_wb(): P must be a square matrix with dimension ncol(Q2) x ncol(Q2).")
  }

  gamma.hat <- as.numeric(gamma.hat)
  expected_len <- nq + np * J
  if (length(gamma.hat) != expected_len) {
    stop(
      paste0(
        "SQBiT_wb(): gamma.hat has length ", length(gamma.hat),
        ", but expected ", expected_len,
        " = ncol(C) + ncol(X) * ncol(Q2)."
      )
    )
  }
  gamma.hat <- matrix(gamma.hat, ncol = 1)

  # ----------------------------------------------------------
  # Sanitize lambda: force a numeric vector of length np
  # ----------------------------------------------------------
  .sqbit_clean_lambda <- function(lambda, np) {
    while (is.list(lambda) && length(lambda) == 1L) {
      lambda <- lambda[[1]]
    }

    if (is.list(lambda) && !is.null(names(lambda))) {
      for (nm in c("lambda.sel", "lambda.star", "lambda.opt", "lambda.min",
                   "lambdac", "lambda", "value")) {
        if (!is.null(lambda[[nm]])) {
          return(.sqbit_clean_lambda(lambda[[nm]], np))
        }
      }
      lambda <- unlist(lambda, recursive = TRUE, use.names = FALSE)
    }

    if (length(lambda) == 1L && is.null(dim(lambda))) {
      return(rep(as.numeric(lambda), np))
    }

    if (is.atomic(lambda) && is.null(dim(lambda))) {
      if (length(lambda) == np) {
        return(as.numeric(lambda))
      }
      if (length(lambda) %% np == 0L) {
        lam.mat <- matrix(as.numeric(lambda), ncol = np, byrow = TRUE)
        return(as.numeric(lam.mat[1, ]))
      }
    }

    lam.mat <- as.matrix(lambda)

    if (ncol(lam.mat) == np) {
      return(as.numeric(lam.mat[1, ]))
    }

    if (nrow(lam.mat) == np && ncol(lam.mat) >= 1L) {
      return(as.numeric(lam.mat[, 1]))
    }

    if (length(lam.mat) %% np == 0L) {
      lam.mat2 <- matrix(as.numeric(lam.mat), ncol = np, byrow = TRUE)
      return(as.numeric(lam.mat2[1, ]))
    }

    stop(
      paste0(
        "SQBiT_wb(): invalid lambda object. Expected a scalar or a vector of length ",
        np, ". Got dim = ", paste(dim(lam.mat), collapse = " x "),
        " and length = ", length(lam.mat), "."
      )
    )
  }

  lambda.vec <- .sqbit_clean_lambda(lambda, np)

  if (!is.numeric(lambda.vec) || length(lambda.vec) != np) {
    stop(
      paste0(
        "SQBiT_wb(): cleaned lambda must be a numeric vector of length ",
        np, ". Got length = ", length(lambda.vec), "."
      )
    )
  }

  # ----------------------------------------------------------
  # Build block penalty directly
  # ----------------------------------------------------------
  kp.blocks <- lapply(seq_len(np), function(j) lambda.vec[j] * P)
  kp0 <- as.matrix(Matrix::bdiag(kp.blocks))

  if (nrow(kp0) != np * J || ncol(kp0) != np * J) {
    stop(
      paste0(
        "SQBiT_wb(): block penalty has dimension ",
        nrow(kp0), " x ", ncol(kp0),
        ", but expected ", np * J, " x ", np * J, "."
      )
    )
  }

  # Construct design matrices
  BQ2 <- eigenMapMatMult(B, Q2)
  W <- as.matrix(kr(X, BQ2, byrow = TRUE))
  W <- cbind(C, W)
  cp.W <- eigenMapMatMult(t(W), W)

  psi.eq <- function(r, eta.j, kernel = kernel) {
    if (kernel == "unif") {
      psi1 <- (eta.j * r - (tau - 1 / 2)) / (eta.j + 1 / (2 * h))
      psi2 <- -(tau - 1) / eta.j + r
      psi3 <- -tau / eta.j + r

      psi1.r <- psi1 * ifelse(abs(psi1) <= h, 1, 0)
      psi2.r <- psi2 * ifelse(psi2 < -h, 1, 0)
      psi3.r <- psi3 * ifelse(psi3 > h, 1, 0)

      psi1.r + psi2.r + psi3.r
    } else if (kernel == "norm") {
      sapply(
        r,
        function(ri) nleqslv::nleqslv(
          0,
          function(u) tau - stats::pnorm(-u / h) + eta.j * (u - ri)
        )$x
      )
    } else {
      stop("SQBiT_wb(): kernel must be 'unif' or 'norm'.")
    }
  }

  zm1 <- matrix(0, nrow = nrow(kp0), ncol = nq)
  zm2 <- matrix(0, nrow = nq, ncol = nrow(kp0) + nq)
  kp <- rbind(zm2, cbind(zm1, kp0))

  if (!all(dim(kp) == dim(cp.W))) {
    stop(
      paste0(
        "SQBiT_wb(): dimension mismatch before chol(). ",
        "dim(kp) = ", paste(dim(kp), collapse = " x "),
        ", dim(crossprod(W)) = ", paste(dim(cp.W), collapse = " x "), "."
      )
    )
  }

  cp.inv <- Matrix::chol2inv(chol((2 / eta.j) * kp + cp.W))
  mat <- eigenMapMatMult(cp.inv, t(W))

  resid <- as.numeric(Y - eigenMapMatMult(W, gamma.hat))

  J.h <- eigenMapMatMult(t(W), W) / n
  f0 <- akj(resid, 0)$dens
  pgrad <- t(
    eigenMapMatMult(
      t(matrix(tau - stats::punif(-resid / h, -1, 1))),
      W
    ) / n
  )

  if (!biascorr) {
    bias <- 0
  } else {
    bias <- as.numeric(
      eigenMapMatMult(W, eigenMapMatMult(AInv(J.h) * (1 / f0), pgrad))
    )
  }

  eta.B <- NULL
  gamma.B <- NULL

  for (b in seq_len(Br)) {
    xis <- rG(n = n, tau = tau)
    Y.b <- eigenMapMatMult(W, gamma.hat) + xis * abs(resid + bias)

    new.u <- matrix((1 / eta.j) * as.numeric(Y.b), ncol = 1)
    new.gamma <- gamma.hat
    new.fit <- eigenMapMatMult(W, new.gamma)

    iter.admm <- 0
    repeat {
      iter.admm <- iter.admm + 1

      old.u <- new.u
      old.fit <- new.fit

      r <- as.numeric(Y.b - old.fit - old.u)
      new.psi <- matrix(psi.eq(r, kernel = kernel, eta.j = eta.j), ncol = 1)

      new.gamma <- eigenMapMatMult(mat, Y.b - new.psi - old.u)

      new.fit <- eigenMapMatMult(W, new.gamma)
      res.new <- Y.b - new.fit
      new.u <- old.u + (new.psi - res.new)

      r.prime <- res.new - new.psi
      r.dual  <- eta.j * (new.fit - old.fit)

      eps.prime <- sqrt(length(r.prime)) * eps.abs +
        eps.rel * max(
          norm(as.matrix(new.psi), type = "2"),
          norm(as.matrix(new.fit), type = "2"),
          norm(as.matrix(Y), type = "2")
        )

      eps.dual <- sqrt(length(r.dual)) * eps.abs +
        eps.rel * norm(as.matrix(new.u), type = "2")

      n.prime <- norm(as.matrix(r.prime), type = "2")
      n.dual  <- norm(as.matrix(r.dual), type = "2")

      if (var.j) {
        old.eta.j <- eta.j

        if (n.prime > zeta * n.dual) {
          eta.j <- incr * eta.j
        } else if (zeta * n.prime < n.dual) {
          eta.j <- eta.j / incr
        }

        if (eta.j != old.eta.j) {
          cp.inv <- Matrix::chol2inv(chol((2 / eta.j) * kp + cp.W))
          mat <- eigenMapMatMult(cp.inv, t(W))
        }
      }

      if ((n.prime < eps.prime && n.dual < eps.dual) ||
          iter.admm >= max.iter) {
        break
      }
    }

    eta.B <- rbind(eta.B, as.numeric(new.gamma[seq_len(nq), , drop = FALSE]))
    gamma.B <- rbind(gamma.B, as.numeric(new.gamma[-seq_len(nq), , drop = FALSE]))
  }

  eta.hat <- as.numeric(gamma.hat[seq_len(nq), , drop = FALSE])
  cis2 <- NULL

  for (j in seq_len(nq)) {
    ci <- stats::quantile(eta.B[, j], c(level / 2, 1 - level / 2))
    cis2 <- rbind(cis2, c(2 * eta.hat[j] - ci[2], 2 * eta.hat[j] - ci[1]))
  }

  if (compute.vc) {
    gamma.hat.vc <- gamma.hat[-seq_len(nq), , drop = FALSE]

    betas.lb <- NULL
    betas.ub <- NULL
    theta.lb <- NULL
    theta.ub <- NULL
    theta.B  <- NULL

    for (j in seq_len(np)) {
      idx <- ((j - 1) * J + 1):(j * J)

      beta.hat.j <- BQ2.eva %*% gamma.hat.vc[idx, , drop = FALSE]
      betas.B.j  <- BQ2.eva %*% t(gamma.B[, idx, drop = FALSE])

      theta.j   <- Q2 %*% gamma.hat.vc[idx, , drop = FALSE]
      theta.B.j <- Q2 %*% t(gamma.B[, idx, drop = FALSE])
      theta.B   <- rbind(theta.B, theta.B.j)

      lb.j <- apply(betas.B.j, 1, stats::quantile, probs = level / 2)
      ub.j <- apply(betas.B.j, 1, stats::quantile, probs = 1 - level / 2)

      theta.lb.j <- apply(theta.B.j, 1, stats::quantile, probs = level / 2)
      theta.ub.j <- apply(theta.B.j, 1, stats::quantile, probs = 1 - level / 2)

      betas.lb <- cbind(betas.lb, 2 * beta.hat.j - matrix(ub.j, ncol = 1))
      betas.ub <- cbind(betas.ub, 2 * beta.hat.j - matrix(lb.j, ncol = 1))

      theta.lb <- cbind(theta.lb, 2 * theta.j - matrix(theta.ub.j, ncol = 1))
      theta.ub <- cbind(theta.ub, 2 * theta.j - matrix(theta.lb.j, ncol = 1))
    }
  } else {
    theta.B  <- NULL
    theta.lb <- NULL
    theta.ub <- NULL
    betas.lb <- NULL
    betas.ub <- NULL
  }

  end_time <- Sys.time()

  list(
    etas = eta.B,
    cis2 = cis2,
    theta.B = theta.B,
    theta.lb = theta.lb,
    theta.ub = theta.ub,
    betas.lb = betas.lb,
    betas.ub = betas.ub,
    result_time = end_time - start_time
  )
}
