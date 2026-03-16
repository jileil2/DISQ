#' Distributed GACV evaluation for a DISQ model
#'
#' @description
#' `SQBiT_int_gacv()` computes the **distributed generalized approximate
#' cross-validation (GACV)** criterion for a given DISQ model specification.
#'
#' The function evaluates the model fit by applying DISQ on each spatial
#' subregion (local machine), computing local residuals and effective
#' degrees of freedom, and aggregating them to form the distributed GACV.
#' This criterion is used for model selection, for example in
#' `forward.disq()`, to determine which covariates should have
#' spatially varying coefficients.
#'
#' Let \eqn{m} denote the number of subregions (or local machines).
#' For observations \eqn{i = 1,\ldots,n_k} on machine \eqn{k},
#' the local residuals are
#'
#' \deqn{
#' \hat{\varepsilon}_{k,i}
#' =
#' Y_{k,i}
#' -
#' \boldsymbol{C}_{k,i}^{\top}\hat{\boldsymbol{\eta}}_{k}
#' -
#' \sum_{j=0}^p
#' X_{k,i,j}\hat{\beta}_{k,j}(\boldsymbol{S}_{k,i}).
#' }
#'
#' The distributed GACV is then computed as
#'
#' \deqn{
#' \mathrm{GACV}
#' =
#' \frac{1}{N - D}
#' \sum_{k=1}^{m}
#' \sum_{i=1}^{n_k}
#' \rho_{\tau}(\hat{\varepsilon}_{k,i}),
#' }
#'
#' where \eqn{D} is the total effective degrees of freedom aggregated
#' across all machines.
#'
#' @param tau Quantile level.
#' @param m Integer. Number of spatial subregions (or local machines).
#' @param maxTt Integer. Reserved for compatibility with other routines.
#' @param Y Numeric response vector of length \eqn{N}.
#' @param C Numeric matrix of covariates with constant coefficients.
#' @param X Numeric matrix of covariates with varying coefficients.
#' @param Q2 Basis transformation matrix used in the spline representation.
#' @param max.iter Maximum number of iterations used in the local fitting step.
#' @param cs Numeric vector controlling local bandwidth scaling for each
#'   subregion.
#' @param eta.j Numeric tuning parameter passed to `SQBiT::SQBiT_gacv()`.
#' @param d Integer. Degree of spline basis.
#' @param r Integer. Smoothness order of the spline penalty.
#' @param Tr Integer matrix with 3 columns defining the triangulation
#'   connectivity.
#' @param V Numeric matrix of vertex coordinates.
#' @param S Numeric matrix of spatial locations for observations.
#' @param triangles Optional list defining the domain decomposition.
#'   Each element contains triangle indices forming one subregion.
#'   If `FALSE`, a simple automatic partition of the triangulation
#'   into `m` groups is constructed.
#' @param nlambda Integer. Number of candidate smoothing parameters used
#'   in the local GACV tuning step.
#' @param new.nlambda Integer. Number of refinement steps in lambda tuning.
#' @param lambda.scale Numeric scaling factor used when expanding
#'   the lambda search grid.
#'
#' @details
#' The function operates in three stages.
#'
#' **Step 1. Domain decomposition**
#'
#' If the argument `triangles` is not supplied, the triangulation is
#' automatically partitioned into `m` groups of triangles. Each group
#' corresponds to one local machine.
#'
#' **Step 2. Local estimation**
#'
#' For each subregion, observations whose spatial locations fall inside
#' the sub-triangulation are collected and used to fit the local model
#' using spline basis functions. The smoothing parameter is selected
#' by minimizing a local GACV criterion through `SQBiT::SQBiT_gacv()`.
#'
#' **Step 3. Aggregation**
#'
#' Residuals and effective degrees of freedom from all machines are
#' aggregated to compute the distributed GACV value.
#'
#' This function is primarily used as an internal utility for model
#' selection procedures such as `forward.disq()`.
#'
#' @return
#' Numeric scalar representing the distributed GACV value.
#'
#' @examples
#' \dontrun{
#' gacv <- SQBiT_int_gacv(
#'   tau = 0.5,
#'   m = length(triangles),
#'   Y = Y,
#'   C = C,
#'   X = X,
#'   Q2 = Q2,
#'   cs = rep(7, length(triangles)),
#'   Tr = Tr,
#'   V = V,
#'   S = S,
#'   triangles = triangles
#' )
#' }
#'
#' @export
SQBiT_int_gacv <- function(tau,
                           Y, C, X, Q2,
                           cs, eta.j = 1,
                           d = 3, r = 1,
                           Tr, V, S,
                           triangles,
                           max.iter = 50,
                           nlambda = 10,
                           new.nlambda = 10,
                           lambda.scale = 100) {

  N <- length(Y)
  nt <- length(triangles)

  resid <- vector("list", nt)
  dfvs  <- numeric(nt)

  for (a in seq_len(nt)) {

    ind.a <- inVT(
      V0 = V,
      Tr0 = Tr[triangles[[a]], , drop = FALSE],
      xx = S[, 1],
      yy = S[, 2]
    )$ind.inside

    Y.a <- Y[ind.a]
    X.a <- as.matrix(X[ind.a, , drop = FALSE])
    S.a <- S[ind.a, , drop = FALSE]

    if (is.null(C)) {
      C.a <- NULL
    } else {
      C.a <- as.matrix(C[ind.a, , drop = FALSE])
      if (!is.null(ncol(C.a)) && ncol(C.a) == 0) C.a <- NULL
    }

    n.a <- length(Y.a)
    if (n.a == 0) next

    B0.sub <- basis(
      V = V,
      Tr = Tr[triangles[[a]], , drop = FALSE],
      d = d,
      r = r,
      Z = S.a
    )

    Q2.sub <- as.matrix(B0.sub$Q2)
    B.sub  <- as.matrix(B0.sub$B)
    K.sub  <- as.matrix(B0.sub$K)

    BQ2.sub <- B.sub %*% Q2.sub
    P.sub   <- t(Q2.sub) %*% K.sub %*% Q2.sub

    W.vary <- kr(X.a, BQ2.sub, byrow = TRUE)

    qC <- if (is.null(C.a)) 0 else ncol(C.a)
    h <- cs[a] * tau * (1 - tau) *
      ((qC + ncol(Q2.sub) * ncol(X.a) + log(n.a)) / n.a)^(2 / 5)

    fit.a <- SQBiT::SQBiT_gacv(
      tau = tau,
      Y = Y.a,
      C = C.a,
      X = X.a,
      B = B.sub,
      Q2 = Q2.sub,
      P = P.sub,
      h = h,
      max.iter = max.iter,
      eta.j1 = eta.j,
      eta.j2 = eta.j,
      nlambda = nlambda,
      new.nlambda = new.nlambda,
      lambda.scale = lambda.scale
    )

    nq <- if (is.null(C.a)) 0 else ncol(C.a)

    W.full <- W.vary
    if (!is.null(C.a)) {
      W.full <- cbind(C.a, W.full)
    }

    cp.W <- eigenMapMatMult(t(W.full), W.full)

    lambdac <- fit.a$lambdac
    if (length(lambdac) == 1) {
      kp <- lambdac * P.sub
    } else {
      kp <- cppkp(diag(lambdac), P.sub)
    }

    if (nq > 0) {
      zm1 <- matrix(0, nrow = nrow(kp), ncol = nq)
      zm2 <- matrix(0, nrow = nq, ncol = nrow(kp) + nq)
      kp <- rbind(zm2, cbind(zm1, kp))
    }

    cp.inv <- Matrix::chol2inv(chol((2 / eta.j) * kp + cp.W))
    dfvs[a] <- sum_cpp(diag(eigenMapMatMult(cp.inv, cp.W)))

    const.part <- 0
    if (!is.null(C.a) && !is.null(fit.a$eta)) {
      const.part <- eigenMapMatMult(C.a, matrix(fit.a$eta))
    }

    vary.part <- eigenMapMatMult(
      W.vary,
      matrix(fit.a$gamma, ncol = 1)
    )

    resid[[a]] <- as.numeric(Y.a - const.part - vary.part)
  }

  loss <- rhotau(u = as.vector(unlist(resid)), tau = tau)
  gacv <- sum(loss) / (N - sum(dfvs))

  gacv
}
