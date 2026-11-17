#' Fit local SQBiT models on worker subregions
#'
#' Fits penalized smoothed quantile regressions on the supplied triangle groups
#' and computes local decorrelated constant-coefficient estimates. Raw spline
#' coefficients are assembled in global triangle order and projected onto the
#' global smoothness constraints.
#'
#' @param tau Quantile level in \eqn{(0,1)}.
#' @param Y Numeric response vector of length \eqn{N}.
#' @param C Numeric matrix (\eqn{N \times q}) of covariates with constant effects.
#' @param X Numeric matrix (\eqn{N \times p}) of covariates with varying coefficients.
#' @param Q2 Global \eqn{Q_2} matrix from the QR decomposition of the smoothness matrix
#' used by the triangulation spline basis.
#' @param cs Positive numeric vector of bandwidth multipliers; one per subregion after empty
#' subregions are removed. The bandwidth in subregion \eqn{a} is computed as
#' \deqn{
#' h_a = c_a \, \tau(1-\tau) \left\{\frac{q + p\,\mathrm{dim}(Q_{2,a}) + \log(n_a)}{n_a}\right\}^{2/5}.
#' }
#' @param eta.j Step size parameter used in local SQBiT fitting (typically \code{1} when step size is varied).
#' @param d Degree of Bernstein polynomials for bivariate triangulation splines.
#' @param r Smoothness order for bivariate triangulation splines.
#' @param Tr Integer triangle index matrix (\eqn{n_T \times 3}) for the global triangulation.
#' @param V Numeric vertex coordinate matrix (\eqn{n_V \times 2}) for the global triangulation.
#' @param S Numeric location matrix (\eqn{N \times 2}) giving coordinates for each observation.
#' @param triangles A non-empty list of integer vectors. Each element specifies triangle indices
#' (rows of \code{Tr}) defining a sub-triangulation (subregion). Empty subregions are removed internally.
#' @param proj.H Optional projection matrix used to smooth/project spline coefficients.
#' If \code{NULL}, it is computed internally from \code{smoothness(V, Tr, d, r)}.
#' @param max.iter Maximum number of ADMM iterations used in local SQBiT and decorrelated-score steps.
#' @param nlambda Number of candidate \eqn{\lambda} values in the local GACV search (SQBiT internal).
#' @param new.nlambda Number of refined \eqn{\lambda} values in the local search (SQBiT internal).
#' @param lambda.scale Scale factor controlling the \eqn{\lambda} grid (SQBiT internal).
#' @param ncores Number of parallel workers. If \code{<=1} or only one subregion, runs sequentially.
#' @param verbose Logical; reserved for future logging (currently minimal).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{eta}}{Stacked local constant-effect estimates (\code{fit.a$eta}) across subregions.}
#'   \item{\code{eta.star}}{Stacked decorrelated-score constant-effect estimates across subregions
#'   (from \code{SQBiT_ds}).}
#'   \item{\code{Ys, Xs, Ss, Cs, Ws}}{Lists of per-subregion response, covariates, locations, and design matrices.}
#'   \item{\code{ns}}{Integer vector of subregion sample sizes.}
#'   \item{\code{hs}}{Numeric vector of subregion bandwidths used for convolution smoothing.}
#'   \item{\code{dQ2Xs}}{Numeric vector; per-subregion dimension \eqn{\mathrm{dim}(Q_{2,a}) \times p}.}
#'   \item{\code{gammas}}{List of per-subregion spline coefficient vectors (\code{fit.a$gamma}).}
#'   \item{\code{Gs}}{List of per-subregion kernel-weighted projection matrices from \code{SQBiT_ds}.}
#'   \item{\code{thetas}}{Global stacked basis coefficients by triangle blocks.}
#'   \item{\code{thetas.sm}}{Smoothed/projected \code{thetas} via \code{proj.H}.}
#'   \item{\code{loss}}{Vector of check-loss contributions \code{rhotau(resid, tau)} across all subregions.}
#'   \item{\code{gacv}}{Global GACV-like criterion computed using finite \code{dfvs} only.}
#'   \item{\code{lambdas}}{List of selected \eqn{\lambda} (scalar) per subregion from local GACV search.}
#'   \item{\code{resid}}{List of per-subregion residual vectors \code{Y - C eta - W gamma}.}
#'   \item{\code{vcs}}{Numeric vector of stacked varying-coefficient contributions \code{W gamma}.}
#'   \item{\code{dfvs}}{Numeric vector of per-subregion degrees-of-freedom for varying effects; may contain \code{NA}.}
#'   \item{\code{time}}{Elapsed runtime.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' ## Step 1. Create triangulation
#' boundaries <- matrix(c(0, 0, 0, 1, 1, 1, 1, 0), ncol = 2, byrow = TRUE)
#' tri <- Triangulation::TriMesh(boundaries, n = 13)
#' Tr <- tri$Tr
#' V <- round(tri$V, 3)
#'
#' ## Step 2. Population locations and keep those inside triangulation
#' S_pop <- cbind(runif(10000), runif(10000))
#' ind_inside <- inVT(V0 = V, Tr0 = Tr, xx = S_pop[, 1], yy = S_pop[, 2])$ind.inside
#' S_pop <- S_pop[ind_inside, , drop = FALSE]
#'
#' ## Step 3. Global triangulation spline basis (population)
#' d <- 3
#' r <- 1
#' B0_pop <- BPST::basis(V = V, Tr = Tr, d = d, r = r, Z = as.matrix(S_pop))
#' Q2 <- B0_pop$Q2
#'
#' ## Step 4. Domain decomposition (sub-triangulations)
#' triangles <- DISQ::subregion_const(Tr = Tr, V = V, n.layer = 1)
#' cs <- rep(7, length(triangles))
#'
#' ## Step 5. Simulate observations
#' N <- 20000
#' idx <- sample.int(nrow(S_pop), N, replace = TRUE)
#' S <- S_pop[idx, , drop = FALSE]
#' C <- matrix(runif(N * 3, -1, 1), ncol = 3)
#' X <- matrix(runif(N * 1, -1, 1), ncol = 1)
#'
#' beta0_fun <- function(s1, s2) sin(pi * s1 * s2)
#' beta1_fun <- function(s1, s2) (1 - (1 - 2 * s1)^2) * (1 - (1 - 2 * s2)^2)
#'
#' beta0 <- beta0_fun(S[, 1], S[, 2])
#' beta1 <- beta1_fun(S[, 1], S[, 2])
#' eta <- c(1, 1, 1)
#' Y <- beta0 + X[, 1] * beta1 + as.numeric(C %*% eta) + rnorm(N)
#'
#' ## Step 6. Fit SQBiT_int
#' fit_int <- SQBiT_int(
#'   tau = 0.5,
#'   Y = Y, C = C, X = X, Q2 = Q2,
#'   cs = cs, eta.j = 1, d = d, r = r,
#'   Tr = Tr, V = V, S = S,
#'   triangles = triangles,
#'   ncores = 4
#' )
#'
#' str(fit_int, max.level = 1)
#' }
#'
#' @export
SQBiT_int <- function(
  tau,
  Y, C, X, Q2,
  cs,
  eta.j = 1, d = 3, r = 1,
  Tr, V, S,
  triangles,
  proj.H = NULL,
  max.iter = 50,
  nlambda = 10, new.nlambda = 10, lambda.scale = 100,
  ncores = max(1L, parallel::detectCores() - 1L),
  verbose = FALSE
) {
  start_time <- Sys.time()

  # Required packages (master)
  requireNamespace("SQBiT", quietly = TRUE)
  requireNamespace("DISQ", quietly = TRUE)
  requireNamespace("Matrix", quietly = TRUE)
  requireNamespace("mgcv", quietly = TRUE)
  requireNamespace("BPST", quietly = TRUE)
  requireNamespace("Triangulation", quietly = TRUE)

  # Basic coercions
  Y <- as.numeric(Y)
  N <- length(Y)
  C <- as.matrix(C)
  X <- as.matrix(X)
  S <- as.matrix(S)

  if (!(tau > 0 && tau < 1)) stop("tau must be in (0,1).")
  if (!is.matrix(Tr) || ncol(Tr) != 3) stop("Tr must be an integer matrix with 3 columns.")
  if (!is.matrix(V) || ncol(V) != 2) stop("V must be a numeric matrix with 2 columns.")
  if (!is.matrix(S) || ncol(S) != 2) stop("S must be a numeric matrix with 2 columns (x,y).")
  if (nrow(C) != N) stop("nrow(C) must equal length(Y).")
  if (nrow(X) != N) stop("nrow(X) must equal length(Y).")
  if (nrow(S) != N) stop("nrow(S) must equal length(Y).")

  if (!is.list(triangles) || length(triangles) == 0L) {
    stop("triangles must be a non-empty list of triangle-index vectors.")
  }

  if (!exists("inVT", mode = "function")) stop("Function 'inVT()' not found.")
  if (!exists("basis", mode = "function")) stop("Function 'basis()' not found.")
  if (!exists("eigenMapMatMult", mode = "function")) stop("Function 'eigenMapMatMult()' not found.")
  if (!exists("smoothness", mode = "function")) stop("Function 'smoothness()' not found.")
  if (!exists("armaInv", mode = "function")) stop("Function 'armaInv()' not found.")
  if (!exists("SQBiT_ds", mode = "function")) stop("Function 'SQBiT_ds()' not found.")

  # Remove empty subregions
  nt0 <- length(triangles)
  keep <- logical(nt0)
  for (a in seq_len(nt0)) {
    ind.a <- inVT(
      V0 = V, Tr0 = Tr[triangles[[a]], , drop = FALSE],
      xx = S[, 1], yy = S[, 2]
    )
    keep[a] <- length(ind.a$ind.inside) > 0L
  }
  triangles <- triangles[keep]
  nt <- length(triangles)
  if (nt == 0L) stop("All subregions are empty: no observations fell inside any sub-triangulation.")

  nq_basis <- (d + 1) * (d + 2) / 2

  # Projection matrix
  if (is.null(proj.H)) {
    H <- as.matrix(smoothness(V, Tr, d, r))
    HH <- eigenMapMatMult(H, t(H))
    nH <- nrow(HH)
    bb <- Matrix::chol2inv(chol(HH + 1e-12 * diag(nH)))
    proj.H <- eigenMapMatMult(bb, H)
    proj.H <- eigenMapMatMult(t(H), proj.H)
  } else {
    proj.H <- as.matrix(proj.H)
  }

  if (length(cs) < nt) stop("Length(cs) must be >= number of non-empty subregions.")
  cs <- cs[seq_len(nt)]

  # One subregion
  process_one <- function(a) {
    ind.a <- inVT(
      V0 = V,
      Tr0 = Tr[triangles[[a]], , drop = FALSE],
      xx = S[, 1],
      yy = S[, 2]
    )$ind.inside

    Y.a <- as.numeric(Y[ind.a])
    X.a <- as.matrix(X[ind.a, , drop = FALSE])
    S.a <- as.matrix(S[ind.a, , drop = FALSE])
    C.a <- as.matrix(C[ind.a, , drop = FALSE])
    n.a <- length(Y.a)

    B0.sub <- basis(
      V = V, Tr = Tr[triangles[[a]], , drop = FALSE],
      d = d, r = r, Z = S.a
    )
    Q2.sub <- as.matrix(B0.sub$Q2)
    B.sub <- as.matrix(B0.sub$B)
    K.sub <- as.matrix(B0.sub$K)

    P.sub <- t(Q2.sub) %*% K.sub %*% Q2.sub
    BQ2.sub <- B.sub %*% Q2.sub
    W.a <- as.matrix(kr(X.a, BQ2.sub, byrow = TRUE))

    h <- cs[a] * tau * (1 - tau) *
      ((ncol(C.a) + ncol(Q2.sub) * ncol(X.a) + log(n.a)) / n.a)^(2 / 5)

    fit.a <- SQBiT::SQBiT_gacv(
      tau = tau, Y = Y.a, C = C.a, X = X.a,
      B = B.sub, Q2 = Q2.sub, P = P.sub,
      h = h, max.iter = max.iter,
      eta.j1 = eta.j, eta.j2 = eta.j,
      nlambda = nlambda, new.nlambda = new.nlambda, lambda.scale = lambda.scale
    )

    # Effective degrees of freedom for the local fit.
    pC <- ncol(C.a)
    df_fit <- fit.a$df
    if (is.null(df_fit) || length(df_fit) == 0L || !is.finite(as.numeric(df_fit)[1])) {
      dfv.a <- NA_real_
    } else {
      dfv.a <- as.numeric(df_fit)[1] - pC
    }

    # ---- DS (original) ----
    lambdac <- as.numeric(fit.a$lambdac)[1]

    fit.dcs <- SQBiT_ds(
      tau = tau, Y = Y.a, C = C.a, X = X.a,
      B = B.sub, Q2 = Q2.sub,
      h = h, lambda = lambdac,
      eta.0 = fit.a$eta,
      gamma.0 = fit.a$gamma,
      zeta = 10, max.iter = 200, s.j = 1, var.j = TRUE
    )

    gamma_v <- matrix(as.numeric(fit.a$gamma), ncol = 1)
    eta_v <- matrix(as.numeric(fit.a$eta), ncol = 1)

    vc_a <- as.numeric(eigenMapMatMult(W.a, gamma_v))
    resid_a <- as.numeric(Y.a - vc_a - as.numeric(eigenMapMatMult(C.a, eta_v)))

    theta_blocks <- vector("list", length(triangles[[a]]))
    for (qq in seq_along(triangles[[a]])) {
      tri_id <- triangles[[a]][qq]
      theta_blocks[[qq]] <- list(
        tri_id = tri_id,
        theta  = fit.a$theta[((qq - 1) * nq_basis + 1):(qq * nq_basis), , drop = FALSE]
      )
    }

    list(
      a = a,
      Y = Y.a, X = X.a, S = S.a, C = C.a, W = W.a,
      n = n.a,
      h = h,
      eta0 = fit.a$eta,
      etads = as.numeric(fit.dcs$eta),
      gamma = fit.a$gamma,
      dQ2X = ncol(Q2.sub) * ncol(X.a),
      G = fit.dcs$G,
      dfv = dfv.a,
      lambdac = lambdac,
      resid = resid_a,
      vc = vc_a,
      theta_blocks = theta_blocks
    )
  }

  # Parallel
  if (isTRUE(ncores <= 1L) || nt <= 1L) {
    results <- lapply(seq_len(nt), process_one)
  } else {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(1)

    cl <- parallel::makeCluster(ncores, outfile = "")
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

    lp <- .libPaths()
    parallel::clusterCall(cl, function(p) .libPaths(p), lp)

    parallel::clusterEvalQ(cl, {
      requireNamespace("SQBiT", quietly = TRUE)
      requireNamespace("DISQ", quietly = TRUE)
      requireNamespace("Matrix", quietly = TRUE)
      requireNamespace("mgcv", quietly = TRUE)
      requireNamespace("BPST", quietly = TRUE)
      requireNamespace("Triangulation", quietly = TRUE)

      library(Matrix)
      library(BPST)

      NULL
    })

    parallel::clusterExport(
      cl,
      varlist = c(
        "tau", "Y", "C", "X", "S", "V", "Tr", "triangles", "cs",
        "eta.j", "d", "r", "nq_basis",
        "max.iter", "nlambda", "new.nlambda", "lambda.scale",
        "inVT", "basis", "kr", "eigenMapMatMult", "smoothness", "armaInv",
        "SQBiT_ds",
        "process_one"
      ),
      envir = environment()
    )

    results <- parallel::parLapply(cl, X = seq_len(nt), fun = process_one)
  }

  # Combine worker estimates in triangle order.
  Ys <- lapply(results, `[[`, "Y")
  Xs <- lapply(results, `[[`, "X")
  Ss <- lapply(results, `[[`, "S")
  Cs <- lapply(results, `[[`, "C")
  Ws <- lapply(results, `[[`, "W")

  ns <- vapply(results, `[[`, numeric(1), "n")
  hs <- vapply(results, `[[`, numeric(1), "h")
  dfvs <- vapply(results, function(z) {
    val <- z[["dfv"]]
    if (length(val) == 0L) NA_real_ else as.numeric(val)[1]
  }, numeric(1))
  dQ2Xs <- vapply(results, `[[`, numeric(1), "dQ2X")

  eta.0 <- do.call(rbind, lapply(results, `[[`, "eta0"))
  eta.ds <- do.call(rbind, lapply(results, `[[`, "etads"))
  gammas <- lapply(results, `[[`, "gamma")
  Gs <- lapply(results, `[[`, "G")
  lambdas <- lapply(results, `[[`, "lambdac")
  resid <- lapply(results, `[[`, "resid")
  vcs <- unlist(lapply(results, `[[`, "vc"), use.names = FALSE)

  thetas_list <- vector("list", length = nrow(Tr))
  for (res in results) for (blk in res$theta_blocks) thetas_list[[blk$tri_id]] <- blk$theta
  thetas <- do.call(rbind, thetas_list)
  thetas.sm <- thetas - eigenMapMatMult(proj.H, thetas)

  loss <- rhotau(u = as.vector(unlist(resid, use.names = FALSE)), tau = tau)

  # Use finite dfvs only in the denominator
  dfvs_ok <- dfvs[is.finite(dfvs)]
  gacv <- sum(loss) / (N - sum(dfvs_ok) - ncol(C))

  end_time <- Sys.time()

  list(
    eta = eta.0,
    eta.star = eta.ds,
    Ys = Ys, Ws = Ws, Cs = Cs, ns = ns, Xs = Xs,
    dQ2Xs = dQ2Xs, gammas = gammas,
    Ss = Ss, Gs = Gs, thetas = thetas, thetas.sm = thetas.sm,
    loss = loss,
    gacv = gacv,
    lambdas = lambdas,
    resid = resid,
    hs = hs,
    vcs = vcs,
    dfvs = dfvs,
    time = end_time - start_time
  )
}
