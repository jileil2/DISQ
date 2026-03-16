#' Forward selection for varying-coefficient identification in DISQ
#'
#' @description
#' `forward.disq()` performs a forward selection procedure for determining
#' which covariates should be modeled with **spatially varying coefficients**
#' and which should remain as **constant coefficients** in the DISQ framework.
#'
#' The procedure is based on the distributed generalized approximate
#' cross-validation (GACV) criterion computed by `SQBiT_int_gacv()`.
#' It begins with a baseline model in which all candidate covariates are
#' treated as constant effects, while the intercept is always modeled as a
#' varying coefficient. At each step, the algorithm moves one covariate from
#' the constant set to the varying set, choosing the covariate that yields
#' the largest reduction in GACV. The procedure stops when no further
#' improvement is achieved.
#'
#' Let \eqn{\mathcal{V}} and \eqn{\mathcal{C}} denote the index sets of
#' covariates with varying and constant coefficients, respectively. The
#' baseline model is
#' \deqn{\mathcal{V}^{(0)} = \{0\}, \qquad \mathcal{C}^{(0)} = \{1,\ldots,p\},}
#' where index \eqn{0} corresponds to the intercept. At iteration \eqn{t},
#' for each \eqn{j \in \mathcal{C}^{(t)}}, the algorithm fits the candidate
#' model obtained by moving covariate \eqn{j} from the constant set to the
#' varying set, evaluates its distributed GACV, and selects the covariate
#' producing the smallest criterion value. The update is accepted only if the
#' GACV decreases.
#'
#' @param tau Quantile level. Must be between 0 and 1.
#' @param maxTt Integer. Maximum number of forward-selection steps.
#' @param Y Numeric response vector.
#' @param C Numeric matrix of candidate covariates initially treated as
#'   constant effects. Each column represents one candidate covariate.
#' @param Q2 Basis transformation matrix used in the spline representation.
#' @param max.iter Integer. Maximum number of iterations used in the local
#'   GACV fitting routine.
#' @param cs Numeric vector controlling local bandwidth scaling for each
#'   subregion. Passed to `SQBiT_int_gacv()`.
#' @param eta.j Positive numeric tuning parameter used in the local
#'   GACV fitting routine.
#' @param d Integer. Degree of spline basis.
#' @param r Integer. Smoothness order of the spline penalty.
#' @param Tr Integer matrix with 3 columns. Each row contains the vertex
#'   indices of one triangle in the triangulation.
#' @param V Numeric matrix with 2 columns giving the vertex coordinates.
#' @param S Numeric matrix of spatial locations, with one row per observation.
#' @param triangles List describing the domain decomposition. Each element
#'   contains the triangle indices for one subregion.
#'
#' @details
#' The algorithm proceeds as follows:
#'
#' **Step 1. Baseline model.**
#'
#' Fit the initial model in which all candidate covariates are treated as
#' constant effects and only the intercept is modeled with a varying
#' coefficient.
#'
#' **Step 2. Candidate updates.**
#'
#' For each covariate currently in the constant set, fit a candidate model in
#' which that covariate is moved to the varying set, while all remaining
#' covariates stay in the constant set.
#'
#' **Step 3. GACV comparison.**
#'
#' Compute the distributed GACV for each candidate model and identify the
#' covariate producing the smallest GACV value.
#'
#' **Step 4. Forward update.**
#'
#' If the best candidate reduces the current GACV, move that covariate from
#' the constant set to the varying set and continue. Otherwise, terminate the
#' procedure.
#'
#' The function returns the final partition of covariates into constant and
#' varying components, together with the sequence of GACV values encountered
#' during selection.
#'
#' In the current implementation, the intercept is always included in the
#' varying-coefficient component through `X = matrix(1, nrow = N, ncol = 1)`
#' in the baseline fit and through `X = cbind(1, ...)` in subsequent
#' candidate fits.
#'
#' This function relies on `SQBiT_int_gacv()` for distributed GACV evaluation
#' under the DISQ framework.
#'
#' @return
#' A list with the following components:
#' \describe{
#'   \item{constant}{Integer vector of column indices of `C` selected as
#'   constant effects in the final model.}
#'   \item{varying}{Integer vector of column indices of `C` selected as
#'   varying effects in the final model.}
#'   \item{gacvs}{List containing the baseline GACV value and the candidate
#'   GACV values evaluated at each forward-selection step.}
#'   \item{cri.base}{Final selected GACV value.}
#' }
#'
#' @examples
#' \dontrun{
#' ## ============================================================
#' ## Example:
#' ## Triangulation -> Domain decomposition -> forward.disq
#' ## True model:
#' ##   - 2 varying coefficients total: intercept + X1
#' ##   - 3 constant coefficients: C1, C2, C3
#' ## ============================================================
#'
#' set.seed(123)
#'
#' ## Step 1. Create triangulation
#' boundaries <- matrix(
#'   c(0, 0,
#'     0, 1,
#'     1, 1,
#'     1, 0),
#'   nrow = 4, byrow = TRUE
#' )
#'
#' tri <- TriMesh(boundaries, n = 13)
#' Tr  <- tri$Tr
#' V   <- round(tri$V, 3)
#'
#' ## Step 2. Generate population locations and keep points inside the domain
#' Npop_raw <- 10000
#' S_pop <- cbind(s1 = runif(Npop_raw), s2 = runif(Npop_raw))
#'
#' ind_inside <- inVT(
#'   V0 = V, Tr0 = Tr,
#'   xx = S_pop[, 1], yy = S_pop[, 2]
#' )$ind.inside
#'
#' S_pop_in <- S_pop[ind_inside, , drop = FALSE]
#' Npop <- nrow(S_pop_in)
#'
#' ## Step 3. Define spatially varying coefficients
#' beta0_fun <- function(s1, s2) {
#'   sin(pi * s1 * s2)
#' }
#'
#' beta1_fun <- function(s1, s2) {
#'   (1 - (1 - 2 * s1)^2) * (1 - (1 - 2 * s2)^2)
#' }
#'
#' beta0_pop <- beta0_fun(S_pop_in[, 1], S_pop_in[, 2])
#' beta1_pop <- beta1_fun(S_pop_in[, 1], S_pop_in[, 2])
#'
#' ## Step 4. Build triangulation spline basis
#' d <- 3
#' r <- 1
#'
#' B0_pop <- basis(
#'   V = V,
#'   Tr = Tr,
#'   d = d,
#'   r = r,
#'   Z = as.matrix(S_pop_in)
#' )
#'
#' Q2 <- B0_pop$Q2
#'
#' ## Step 5. Domain decomposition
#' triangles <- subregion_const(Tr = Tr, V = V, n.layer = 1)
#' m <- length(triangles)
#'
#' ## Optional plot
#' sub_tri_plot(triangles, V = V, Tr = Tr, main = expression(Delta))
#'
#' ## Step 6. Simulate observations
#' N <- 20000
#' tau <- 0.5
#'
#' idx <- sample.int(Npop, N, replace = TRUE)
#' S <- S_pop_in[idx, , drop = FALSE]
#'
#' ## Three true constant covariates
#' pC <- 3
#' C.true <- matrix(runif(N * pC, -1, 1), nrow = N, ncol = pC)
#' colnames(C.true) <- paste0("C", 1:pC)
#' eta.true <- c(1.0, -1.0, 0.5)
#'
#' ## Candidate covariates for forward selection
#' ## Only X1 truly has a varying coefficient
#' pX <- 3
#' X <- matrix(runif(N * pX, -1, 1), nrow = N, ncol = pX)
#' colnames(X) <- paste0("X", 1:pX)
#'
#' beta0 <- beta0_pop[idx]
#' beta1 <- beta1_pop[idx]
#'
#' Y <- beta0 +
#'   X[, 1] * beta1 +
#'   as.numeric(C.true %*% eta.true) +
#'   rnorm(N)
#'
#' ## Step 7. Forward selection using distributed GACV
#' ## The intercept is always varying internally.
#' fit.sel <- forward.disq(
#'   tau = tau,
#'   Y = Y,
#'   C = X,
#'   Q2 = Q2,
#'   cs = rep(7, m),
#'   eta.j = 0.32,
#'   d = d,
#'   r = r,
#'   triangles = triangles,
#'   Tr = Tr,
#'   V = V,
#'   S = S
#' )
#'
#' ## Selected candidate covariates
#' fit.sel$constant
#' fit.sel$varying
#' fit.sel$cri.base
#' }
#'
#' @export
forward.disq <- function(tau, maxTt = 100,
                         Y, C, Q2, max.iter = 50,
                         cs, eta.j = 1, d = 3, r = 1,
                         Tr, V, S, triangles) {

  p <- ncol(C)
  N <- length(Y)

  base.gacv <- SQBiT_int_gacv(
    tau = tau,
    Y = Y,
    C = C,
    X = matrix(1, nrow = N, ncol = 1),
    Q2 = Q2,
    cs = cs,
    eta.j = eta.j,
    d = d,
    r = r,
    Tr = Tr,
    V = V,
    S = S,
    triangles = triangles,
    max.iter = max.iter
  )

  cri.base <- base.gacv
  index.C <- seq_len(p)
  index.X <- integer(0)

  gacvs <- vector("list", length = p + 1)
  gacvs[[1]] <- base.gacv
  t <- 1

  while (length(index.C) >= 1 && t < maxTt) {
    cris <- numeric(length(index.C))

    for (jj in seq_along(index.C)) {
      j <- index.C[jj]

      keep.C <- setdiff(index.C, j)

      C.j <- if (length(keep.C) > 0) {
        C[, keep.C, drop = FALSE]
      } else {
        NULL
      }

      X.j <- cbind(1, C[, c(index.X, j), drop = FALSE])

      gacv.j <- SQBiT_int_gacv(
        tau = tau,
        Y = Y,
        C = C.j,
        X = X.j,
        Q2 = Q2,
        cs = cs,
        eta.j = eta.j,
        d = d,
        r = r,
        Tr = Tr,
        V = V,
        S = S,
        triangles = triangles,
        max.iter = max.iter
      )

      cris[jj] <- gacv.j
    }

    t <- t + 1
    gacvs[[t]] <- cris

    jbest.local <- which.min(cris)
    if (cris[jbest.local] > cri.base) break

    jbest <- index.C[jbest.local]
    index.X <- c(index.X, jbest)
    index.C <- index.C[-jbest.local]
    cri.base <- cris[jbest.local]
  }

  list(
    constant = sort(as.numeric(index.C)),
    varying = sort(as.numeric(index.X)),
    gacvs = gacvs[1:t],
    cri.base = cri.base
  )
}
