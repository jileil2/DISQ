# ============================================================
# disq_vc_wb(): distributed wild bootstrap for varying coefficients
# ============================================================

#' Distributed wild-bootstrap inference for varying coefficients in DISQ
#'
#' Performs distributed wild-bootstrap inference for the varying-coefficient
#' component in the DISQ framework. For each nonempty subregion, the function
#' fits a local wild-bootstrap smoothed quantile regression model using
#' \code{SQBiT_wb()}, reconstructs the bootstrap spline coefficients on the
#' global triangulation, optionally applies the global projection smoother,
#' and constructs pointwise confidence intervals for the varying coefficients
#' evaluated on \code{B.eva}.
#'
#' @param tau Quantile level in \eqn{(0,1)}.
#' @param lambdas A list of region-specific smoothing parameters aligned with
#'   the local subregions.
#' @param theta.sm A matrix of globally smoothed spline coefficients for the
#'   varying-coefficient component.
#' @param Ys A list of response vectors for the nonempty subregions.
#' @param Cs A list of design matrices for constant-effect covariates.
#' @param Xs A list of design matrices for varying-effect covariates.
#' @param Ss A list of spatial location matrices for the nonempty subregions.
#' @param etas A matrix of local constant-coefficient estimates aligned with
#'   the supplied subregions. If an intercept is included, it is assumed to be
#'   the first column and is dropped when not needed by \code{C}.
#' @param gammas A list of local spline coefficient vectors aligned with the
#'   supplied subregions.
#' @param B.eva Basis evaluation matrix used to evaluate the reconstructed
#'   varying coefficients.
#' @param Br Number of wild-bootstrap replications.
#' @param max.iter Maximum number of ADMM iterations passed to \code{SQBiT_wb()}.
#' @param cs Region-specific bandwidth scaling constants aligned with the
#'   supplied local subregions.
#' @param eta.j Initial ADMM tuning parameter passed to \code{SQBiT_wb()}.
#' @param d Degree of the spline basis.
#' @param r Smoothness order for the spline basis.
#' @param level Significance level for pointwise confidence intervals.
#' @param Tr Global triangulation matrix.
#' @param V Global vertex matrix.
#' @param triangles A list of triangle-index vectors aligned with the local
#'   nonempty subregions.
#' @param proj.H Projection operator used for global smoothing of reconstructed
#'   bootstrap spline coefficients. If \code{NULL}, the projection step is skipped.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{beta.ub}}{Upper pointwise confidence bounds evaluated on
#'   \code{B.eva}.}
#'   \item{\code{beta.lb}}{Lower pointwise confidence bounds evaluated on
#'   \code{B.eva}.}
#'   \item{\code{beta.est}}{Estimated varying coefficients evaluated on
#'   \code{B.eva}.}
#'   \item{\code{time}}{Elapsed computation time.}
#' }
#'
#' @export
disq_vc_wb <- function(tau, lambdas, theta.sm,
                       Ys, Cs, Xs, Ss, etas, gammas, B.eva,
                       Br = 500, max.iter = 50,
                       cs, eta.j = 1, d = 3, r = 1, level = 0.05,
                       Tr, V, triangles, proj.H = NULL) {
  start_time <- Sys.time()

  # ----------------------------------------------------------
  # helper: extract one local lambda vector of length p
  # ----------------------------------------------------------
  .disq_extract_lambda <- function(obj, p, region_id = NA_integer_) {
    nm_try <- c(
      "lambda.sel", "lambda.star", "lambda.opt", "lambda.min",
      "lambdac", "lambda", "value"
    )

    if (is.list(obj) && !is.data.frame(obj)) {
      if (!is.null(names(obj))) {
        for (nm in nm_try) {
          if (!is.null(obj[[nm]])) {
            return(.disq_extract_lambda(obj[[nm]], p, region_id))
          }
        }
      }
      if (length(obj) == 1L) {
        return(.disq_extract_lambda(obj[[1]], p, region_id))
      }
      obj <- unlist(obj, recursive = TRUE, use.names = FALSE)
    }

    if (length(obj) == 1L) {
      return(matrix(rep(as.numeric(obj), p), nrow = 1))
    }

    if (is.atomic(obj) && is.null(dim(obj)) && length(obj) == p) {
      return(matrix(as.numeric(obj), nrow = 1))
    }

    if (is.atomic(obj) && is.null(dim(obj)) && length(obj) %% p == 0L) {
      mat <- matrix(as.numeric(obj), ncol = p, byrow = TRUE)
      return(matrix(mat[1, ], nrow = 1))
    }

    mat <- as.matrix(obj)

    if (ncol(mat) == p) {
      return(matrix(as.numeric(mat[1, ]), nrow = 1))
    }

    if (nrow(mat) == p && ncol(mat) >= 1L) {
      return(matrix(as.numeric(mat[, 1]), nrow = 1))
    }

    if (length(mat) %% p == 0L) {
      mat2 <- matrix(as.numeric(mat), ncol = p, byrow = TRUE)
      return(matrix(mat2[1, ], nrow = 1))
    }

    stop(
      paste0(
        "disq_vc_wb(): cannot extract a valid local lambda vector in region ",
        region_id, ". Need a scalar or a selected vector of length ", p,
        ". Got object with dim = ",
        paste(dim(mat), collapse = " x "),
        " and length = ", length(mat), "."
      )
    )
  }

  nreg <- length(Ys)
  nq <- (d + 1) * (d + 2) / 2

  if (length(Cs) != nreg || length(Xs) != nreg || length(Ss) != nreg ||
    length(gammas) != nreg || length(lambdas) != nreg ||
    length(triangles) != nreg || length(cs) != nreg) {
    stop(
      "disq_vc_wb(): Ys, Cs, Xs, Ss, gammas, lambdas, triangles, and cs must all have the same length."
    )
  }

  if (!is.matrix(etas)) {
    etas <- as.matrix(etas)
  }
  if (nrow(etas) != nreg) {
    stop("disq_vc_wb(): `etas` must have one row for each local subregion.")
  }

  B.eva <- as.matrix(B.eva)
  theta.sm <- as.matrix(theta.sm)

  thetas.B <- vector("list", Br)
  for (b in seq_len(Br)) {
    thetas.B[[b]] <- vector("list", nrow(Tr))
  }

  # ----------------------------------------------------------
  # Local bootstrap fitting and reconstruction
  # ----------------------------------------------------------
  for (a in seq_len(nreg)) {
    Y.a <- Ys[[a]]
    C.a <- as.matrix(Cs[[a]])
    X.a <- as.matrix(Xs[[a]])
    S.a <- Ss[[a]]
    tri.a <- triangles[[a]]

    if (length(Y.a) == 0) next

    n.a <- length(Y.a)

    B0.sub <- basis(
      V = V,
      Tr = Tr[tri.a, , drop = FALSE],
      d = d,
      r = r,
      Z = S.a
    )

    Q2.sub <- as.matrix(B0.sub$Q2)
    B.sub <- as.matrix(B0.sub$B)
    BQ2.sub <- as.matrix(B.sub %*% Q2.sub)

    K.sub <- B0.sub$K
    P.sub <- t(Q2.sub) %*% K.sub %*% Q2.sub

    W.a <- as.matrix(kr(X.a, BQ2.sub, byrow = TRUE))

    h <- cs[a] * tau * (1 - tau) *
      ((ncol(C.a) + ncol(W.a) + log(n.a)) / n.a)^(2 / 5)

    nq.const <- ncol(C.a)
    eta.a <- as.numeric(etas[a, ])

    if (length(eta.a) == nq.const + 1L) {
      eta.a <- eta.a[-1]
    } else if (length(eta.a) != nq.const) {
      stop(
        paste0(
          "disq_vc_wb(): local constant coefficient length mismatch in region ", a,
          ". length(etas[a, ]) = ", length(as.numeric(etas[a, ])),
          ", but ncol(C.a) = ", nq.const, "."
        )
      )
    }

    gamma.a <- c(eta.a, as.numeric(gammas[[a]]))

    p <- ncol(X.a)
    lambda.a <- .disq_extract_lambda(lambdas[[a]], p = p, region_id = a)

    fit.wb.a <- SQBiT_wb(
      h = h,
      tau = tau,
      Y = Y.a,
      C = C.a,
      X = X.a,
      P = P.sub,
      B = B.sub,
      Q2 = Q2.sub,
      max.iter = max.iter,
      lambda = lambda.a,
      eps.abs = 1e-4,
      eps.rel = 1e-2,
      kernel = "unif",
      var.j = FALSE,
      zeta = 10,
      incr = 2,
      eta.j = eta.j,
      BQ2.eva = BQ2.sub,
      compute.vc = TRUE,
      Br = Br,
      level = level,
      gamma.hat = gamma.a,
      biascorr = TRUE
    )

    theta.Bs <- fit.wb.a$theta.B
    if (is.null(theta.Bs)) {
      stop("disq_vc_wb(): `SQBiT_wb()` returned NULL `theta.B`.")
    }

    for (b in seq_len(Br)) {
      bbb <- matrix(theta.Bs[, b], ncol = p, byrow = FALSE)

      for (qq in seq_along(tri.a)) {
        idx <- ((qq - 1) * nq + 1):(qq * nq)
        thetas.B[[b]][[tri.a[qq]]] <- bbb[idx, , drop = FALSE]
      }
    }
  }

  # ----------------------------------------------------------
  # Global smoothing of bootstrap coefficients
  # ----------------------------------------------------------
  theta.B.sm <- vector("list", Br)
  for (b in seq_len(Br)) {
    if (any(vapply(thetas.B[[b]], is.null, logical(1)))) {
      stop(
        "disq_vc_wb(): some triangle bootstrap blocks are missing. Check alignment between local subregions and `triangles`."
      )
    }

    theta.b <- do.call(rbind, thetas.B[[b]])
    theta.b <- as.matrix(theta.b)

    if (is.null(proj.H)) {
      theta.B.sm[[b]] <- theta.b
    } else {
      proj.H.use <- as.matrix(proj.H)
      theta.B.sm[[b]] <- theta.b - eigenMapMatMult(proj.H.use, theta.b)
    }
  }

  # ----------------------------------------------------------
  # Evaluate point estimates and confidence intervals
  # ----------------------------------------------------------
  beta.est <- NULL
  beta.lb <- NULL
  beta.ub <- NULL

  for (j in seq_len(ncol(theta.sm))) {
    beta.est.B.j <- NULL

    for (b in seq_len(Br)) {
      beta.est.B.j <- cbind(
        beta.est.B.j,
        eigenMapMatMult(B.eva, theta.B.sm[[b]][, j, drop = FALSE])
      )
    }

    beta.est.j <- eigenMapMatMult(B.eva, theta.sm[, j, drop = FALSE])

    beta.est <- cbind(beta.est, beta.est.j)
    beta.ub <- cbind(
      beta.ub,
      2 * beta.est.j - apply(beta.est.B.j, 1, quantile, probs = level / 2)
    )
    beta.lb <- cbind(
      beta.lb,
      2 * beta.est.j - apply(beta.est.B.j, 1, quantile, probs = 1 - level / 2)
    )
  }

  end_time <- Sys.time()

  list(
    beta.ub = beta.ub,
    beta.lb = beta.lb,
    beta.est = beta.est,
    time = end_time - start_time
  )
}
