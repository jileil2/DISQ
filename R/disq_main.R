# ============================================================
# disq_main(): wrapper for the full DISQ procedure
# ============================================================

#' Distributed Inference for Spatial Quantile (DISQ)
#'
#' This function performs distributed inference for quantile spatial model
#' (QSM) over a triangulated spatial domain using DISQ. It returns estimates
#' for the constant coefficients, smoothness-corrected varying coefficients,
#' and confidence intervals for the constant coefficients. Optionally, it also
#' performs distributed wild bootstrap for the varying coefficients.
#'
#' DISQ combines local SQBiT estimation on subregions with a communication-based
#' aggregation algorithm. Instead of fitting one global model on the full data,
#' it divides the spatial domain into subregions, estimates local models in
#' parallel, aggregates the varying coefficients through a projection-based
#' smoothness correction, and then refines the constant coefficients through
#' iterative communication across local machines.
#'
#' \strong{What is Quantile Spatial Model (QSM)?}
#'
#' The proposed Quantile Spatial Model (QSM) is given by:
#' \deqn{
#' Q_{\tau}(Y_i \mid \mathbf{C}_i,\mathbf{X}_i, \mathbf{S}_i) =
#' \mathbf{C}_i^{\top}\boldsymbol{\eta}_{\tau}
#' + \sum_{j=0}^p X_{ij}\beta_{j,\tau}(\mathbf{S}_i),
#' \quad i = 1, \ldots, n.
#' }
#' where \eqn{Q_Y(\tau \mid \cdot) = \inf\{y: P(Y \leq y \mid \cdot) \geq \tau\}}
#' denotes the \eqn{\tau}th conditional quantile of \eqn{Y},
#' \eqn{\boldsymbol{\eta}_{\tau}} captures the constant effects, and
#' \eqn{\beta_{j,\tau}(\cdot)} is the unknown spatially varying coefficient
#' function over the spatial domain.
#' \cr
#'
#' This semi-parametric model captures spatial heterogeneity through the
#' spatially varying coefficients and non-spatial effects through the constant
#' coefficients. The model is flexible and interpretable, and DISQ is designed
#' to make estimation scalable for large spatial datasets.
#'
#' \strong{How to use DISQ?}
#'
#' In general, there are five steps to analyze large spatial data using DISQ.
#'
#' 1. Create a triangulation mesh for your spatial domain. This can be done
#'    using package \code{Triangulation}. It can be downloaded from GitHub using
#'    \code{install_github("funstatpackages/Triangulation")}.
#' 2. Generate the bivariate triangulation spline basis given a polynomial degree
#'    \eqn{d} and a smoothness order \eqn{r}. This can be done using package
#'    \code{BPST}. It can be downloaded from GitHub using
#'    \code{install_github("FIRST-Data-Lab/BPST")}.
#' 3. Construct the subregions for distributed estimation. We recommend users
#'    to use \code{subregion_const()} for domain decomposition.
#' 4. If it is unclear which covariates should receive constant coefficients
#'    and which should receive varying coefficients, we recommend users to first
#'    apply \code{forward.disq()} for model selection.
#' 5. Fit the final QSM using \code{disq_main()}.
#'
#' \strong{Practical recommendation}:
#'
#' For domain decomposition, we recommend using \code{subregion_const()} to
#' construct the subregions from the global triangulation. If users do not know
#' which covariates should enter the model with constant effects and which
#' should enter with spatially varying effects, we recommend applying
#' \code{forward.disq()} before calling \code{disq_main()}.
#'
#' \strong{Communication and optimization}:
#'
#' DISQ first performs local SQBiT estimation on each subregion. The varying
#' coefficients are then aggregated through a projection-based smoothness
#' correction so that the final estimates are globally smooth across subregion
#' boundaries. The constant coefficients are further refined through a
#' communication-based multi-round aggregation procedure implemented in
#' \code{disq()}. The optimization in the communication step is based on
#' gradient-based methods such as Barzilai-Borwein gradient descent or Adam.
#'
#' @param tau quantile level
#' @param B.eva optional global basis evaluation matrix. If supplied, the
#'   smoothness-corrected varying coefficients are evaluated at the rows of
#'   \code{B.eva}. If \code{NULL}, \code{betas = NULL} is returned.
#' @param Y response variable
#' @param C covariates with constant effects
#' @param X covariates with varying coefficients. In most applications,
#'   the intercept should be included in \code{X} if an intercept-varying
#'   component is desired.
#' @param Q2 Q2 matrix from the QR decomposition of the smoothness matrix of
#'   triangles
#' @param eta.j step size parameter passed to \code{SQBiT_int()} for local fitting
#' @param triangles list of sub-triangulations defining the domain decomposition.
#'   We recommend constructing \code{triangles} using
#'   \code{subregion_const(Tr = Tr, V = V, ...)}.
#' @param Tr triangulation matrix
#' @param V vertex matrix of the triangulation
#' @param S spatial locations
#' @param maxTt maximum number of communication rounds in \code{disq()}
#' @param method optimization method used in the communication step
#' @param cs positive numeric vector of local bandwidth multipliers, one per
#'   subregion
#' @param master optional index of the master machine. If \code{NULL}, the
#'   function selects the master machine automatically.
#' @param c positive constant used in the communication step in \code{disq()}
#' @param stepMax maximum step size passed to the inner optimizer in \code{disq()}
#' @param level significance level
#' @param alpha learning rate used by Adam
#' @param beta1 exponential decay rate for the first moment in Adam
#' @param beta2 exponential decay rate for the second moment in Adam
#' @param proj.H optional projection matrix for smoothness correction of varying
#'   coefficients
#' @param d degree of Bernstein basis polynomials
#' @param r smoothness order of bivariate triangulation splines
#' @param max.iter maximum number of iterations in local SQBiT fitting
#' @param nlambda number of candidate tuning parameters in local GACV search
#' @param new.nlambda number of refined tuning parameters in local GACV search
#' @param lambda.scale scale factor controlling the local tuning-parameter grid
#' @param ncores number of parallel workers used in \code{SQBiT_int()}
#' @param verbose logical; passed to \code{SQBiT_int()}
#' @param vc.wb whether to perform distributed wild bootstrap for varying
#'   coefficients (TRUE/FALSE)
#' @param Br number of wild-bootstrap replications
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{time}}{Computational time used.}
#'   \item{\code{eta}}{Estimated constant coefficients after communication-based aggregation.}
#'   \item{\code{eta.full}}{Full coefficient vector returned by \code{disq()}.}
#'   \item{\code{betas}}{Estimated varying coefficients evaluated at \code{B.eva}.}
#'   \item{\code{thetas.sm}}{Smoothness-corrected spline coefficients for the varying effects.}
#'   \item{\code{beta.lb}}{Lower confidence limits for varying coefficients from wild bootstrap.}
#'   \item{\code{beta.ub}}{Upper confidence limits for varying coefficients from wild bootstrap.}
#'   \item{\code{fit.vc.wb}}{Output from \code{disq_vc_wb()} if \code{vc.wb = TRUE}.}
#'   \item{\code{cis}}{Confidence intervals for the constant coefficients.}
#'   \item{\code{ci.trace}}{Trace of confidence intervals across communication rounds.}
#'   \item{\code{eta.trace}}{Trace of constant coefficient estimates across communication rounds.}
#'   \item{\code{fit.int}}{Output from \code{SQBiT_int()}.}
#'   \item{\code{fit.disq}}{Output from \code{disq()}.}
#'   \item{\code{master}}{Index of the selected master machine.}
#'   \item{\code{m}}{Number of nonempty subregions used in fitting.}
#' }
#'
#' @examples
#' \donttest{
#' ###### Step 1. Create a Triangulation Mesh
#' boundaries <- matrix(c(0, 0,
#'                        0, 1,
#'                        1, 1,
#'                        1, 0), nrow = 4, byrow = TRUE)
#' tri <- TriMesh(boundaries, n = 13)
#' Tr <- tri$Tr
#' V <- round(tri$V, 3)
#'
#' # Population locations
#' S_pop <- cbind(s1 = runif(10000), s2 = runif(10000))
#'
#' # Determine which points are inside triangulation
#' ind_inside <- inVT(V0 = V, Tr0 = Tr, xx = S_pop[, 1], yy = S_pop[, 2])$ind.inside
#' S_pop_in <- S_pop[ind_inside, , drop = FALSE]
#' Npop <- nrow(S_pop_in)
#'
#' ###### Step 2. Generate bivariate triangulation basis
#' d <- 3
#' r <- 1
#' B0_pop <- basis(V = V, Tr = Tr, d = d, r = r, Z = as.matrix(S_pop_in))
#' Q2 <- B0_pop$Q2
#' B_pop <- B0_pop$B
#' K <- B0_pop$K
#' P <- t(Q2) %*% K %*% Q2
#'
#' ###### Step 3. Construct subregions using subregion_const()
#' triangles <- subregion_const(Tr = Tr, V = V, n.layer = 1)
#' m <- length(triangles)
#'
#' ###### Step 4. Define coefficient functions and simulate data
#' beta0_fun <- function(s1, s2) sin(pi * s1 * s2)
#' beta1_fun <- function(s1, s2) (1 - (1 - 2 * s1)^2) * (1 - (1 - 2 * s2)^2)
#'
#' beta0_pop <- beta0_fun(S_pop_in[, 1], S_pop_in[, 2])
#' beta1_pop <- beta1_fun(S_pop_in[, 1], S_pop_in[, 2])
#'
#' n <- 20000
#' tau <- 0.5
#' idx <- sample.int(Npop, n, replace = TRUE)
#' S <- S_pop_in[idx, , drop = FALSE]
#'
#' # Constant covariates
#' C0 <- matrix(runif(3 * n, -1, 1), ncol = 3)
#' colnames(C0) <- paste0("C", 1:3)
#' eta.true <- c(1, 1, 1)
#'
#' # Candidate covariates for model selection
#' Xcand <- matrix(runif(3 * n, -1, 1), ncol = 3)
#' colnames(Xcand) <- paste0("X", 1:3)
#'
#' beta0 <- beta0_pop[idx]
#' beta1 <- beta1_pop[idx]
#'
#' Y <- beta0 + Xcand[, 1] * beta1 + as.numeric(C0 %*% eta.true) + rnorm(n)
#'
#' ###### Step 5. Use forward.disq() if the roles of covariates are unclear
#' cs <- rep(7, m)
#'
#' cv.sets <- forward.disq(
#'   tau = tau,
#'   Y = Y,
#'   C = Xcand,
#'   Q2 = Q2,
#'   cs = cs,
#'   eta.j = 0.32,
#'   d = d,
#'   r = r,
#'   triangles = triangles,
#'   Tr = Tr,
#'   V = V,
#'   S = S
#' )
#'
#' constant.sets <- cv.sets$constant
#' constant.sets <- constant.sets[order(constant.sets)]
#' constant.sets <- as.numeric(constant.sets)
#'
#' varying.sets <- cv.sets$varying
#' varying.sets <- varying.sets[order(varying.sets)]
#' varying.sets <- as.numeric(varying.sets)
#'
#' ###### Step 6. Build the final design matrices
#' C.final <- cbind(C0, Xcand[, constant.sets, drop = FALSE])
#' X.final <- Xcand[, varying.sets, drop = FALSE]
#'
#' # Keep the intercept in the varying part
#' if (length(varying.sets) == 0) {
#'   X.final <- matrix(1, nrow = n, ncol = 1)
#' } else {
#'   X.final <- cbind(1, X.final)
#' }
#'
#' ###### Step 7. Evaluate the basis if varying coefficients are to be returned
#' B.eva <- as.matrix(B_pop[idx, ] %*% Q2)
#'
#' ###### Step 8. Fit DISQ
#' fit <- disq_main(
#'   tau = tau,
#'   B.eva = B.eva,
#'   Y = Y,
#'   C = C.final,
#'   X = X.final,
#'   Q2 = Q2,
#'   eta.j = 0.32,
#'   triangles = triangles,
#'   Tr = Tr,
#'   V = V,
#'   S = S,
#'   maxTt = 15,
#'   method = "GD-BB",
#'   cs = cs,
#'   c = 2.5,
#'   stepMax = 1e-3,
#'   level = 0.05,
#'   alpha = 1e-5,
#'   beta1 = 0.9,
#'   beta2 = 0.999,
#'   d = d,
#'   r = r,
#'   max.iter = 50,
#'   nlambda = 10,
#'   new.nlambda = 10,
#'   lambda.scale = 100,
#'   ncores = 4,
#'   verbose = FALSE,
#'   vc.wb = FALSE
#' )
#'
#' fit$eta
#' }
#'
#' @export
#' @references
#' Kim, M., Wang, L., and Wang, H. J. (2025). Estimation and inference of
#' quantile spatially varying coefficient models over complicated domains.
#' \emph{Journal of the American Statistical Association}, 120(551), 1853--1867.
#' Taylor \& Francis.
#'
#' @references
#' Lin, J., Wang, H. J., Wang, L., and Kim, M. (2026). Smoothed quantile
#' regression for spatial data. \emph{Journal of Computational and Graphical
#' Statistics}, just-accepted, 1--21. Taylor \& Francis.
#'
#' @references
#' Boyd, S., Parikh, N., Chu, E., Peleato, B., Eckstein, J., and others
#' (2011). Distributed optimization and statistical learning via the alternating
#' direction method of multipliers. \emph{Foundations and Trends in Machine
#' Learning}, 3(1), 1--122.
disq_main <- function(
    tau = 0.5,
    B.eva = NULL,
    Y, C, X, Q2, eta.j = 0.32,
    triangles,
    Tr, V, S,
    maxTt = 15,
    method = "GD-BB",
    cs = rep(7, length(triangles)),
    master = NULL,
    c = 2.5,
    stepMax = 1e-3,
    level = 0.05,
    alpha = 1e-5,
    beta1 = 0.9,
    beta2 = 0.999,
    proj.H = NULL,
    d = 3,
    r = 1,
    max.iter = 50,
    nlambda = 10,
    new.nlambda = 10,
    lambda.scale = 100,
    ncores = max(1L, parallel::detectCores() - 1L),
    verbose = FALSE,
    vc.wb = FALSE,
    Br = 500
) {
  start_time <- Sys.time()

  C <- as.matrix(C)
  X <- as.matrix(X)
  q <- ncol(C)
  p <- ncol(X)

  if (length(cs) != length(triangles)) {
    stop("disq_main(): `cs` must have the same length as `triangles`.")
  }

  # ------------------------------------------------------------
  # Step 1. Local estimation and internal decomposition
  # ------------------------------------------------------------
  fit_int <- SQBiT_int(
    tau = tau,
    Y = Y,
    C = C,
    X = X,
    Q2 = Q2,
    cs = cs,
    eta.j = eta.j,
    d = d,
    r = r,
    Tr = Tr,
    V = V,
    S = S,
    triangles = triangles,
    proj.H = proj.H,
    max.iter = max.iter,
    nlambda = nlambda,
    new.nlambda = new.nlambda,
    lambda.scale = lambda.scale,
    ncores = ncores,
    verbose = verbose
  )

  Ys     <- fit_int$Ys
  Cs     <- fit_int$Cs
  Xs     <- fit_int$Xs
  Ss     <- fit_int$Ss
  Ws     <- fit_int$Ws
  Gs     <- fit_int$Gs
  gammas <- fit_int$gammas
  eta.0  <- fit_int$eta.star

  if (!is.null(fit_int$lambdac)) {
    lambdas <- fit_int$lambdac
  } else if (!is.null(fit_int$lambda.sel)) {
    lambdas <- fit_int$lambda.sel
  } else if (!is.null(fit_int$lambdas)) {
    lambdas <- fit_int$lambdas
  } else {
    lambdas <- NULL
  }

  m <- length(Ys)

  # ------------------------------------------------------------
  # Step 2. Aggregate varying coefficients through smoothness correction
  # ------------------------------------------------------------
  thetas.sm <- fit_int$thetas.sm

  if (is.null(B.eva)) {
    betas <- NULL
  } else {
    B.eva <- as.matrix(B.eva)
    betas <- vector("list", p)
    for (j in seq_len(p)) {
      betas[[j]] <- as.numeric(
        eigenMapMatMult(B.eva, matrix(thetas.sm[, j], ncol = 1))
      )
    }
  }

  # ------------------------------------------------------------
  # Step 3. Optional projected wild bootstrap for varying coefficients
  # ------------------------------------------------------------
  fit.vc.wb <- NULL
  beta.lb <- NULL
  beta.ub <- NULL

  if (isTRUE(vc.wb)) {
    if (is.null(B.eva)) {
      stop("disq_main(): `B.eva` must be supplied when `vc.wb = TRUE`.")
    }
    if (is.null(lambdas)) {
      stop("disq_main(): local smoothing parameters are unavailable for `vc.wb = TRUE`.")
    }
    if (is.null(gammas)) {
      stop("disq_main(): local fitted spline coefficients `gammas` are unavailable for `vc.wb = TRUE`.")
    }
    if (is.null(eta.0)) {
      stop("disq_main(): local constant coefficients `eta.0` are unavailable for `vc.wb = TRUE`.")
    }

    if (!is.null(fit_int$triangles.nonempty)) {
      triangles.wb <- fit_int$triangles.nonempty
    } else if (!is.null(fit_int$triangles.use)) {
      triangles.wb <- fit_int$triangles.use
    } else if (length(triangles) == length(Ys)) {
      triangles.wb <- triangles
    } else {
      stop(
        paste0(
          "disq_main(): cannot align `triangles` with nonempty local subregions. ",
          "Please have `SQBiT_int()` return `triangles.nonempty` or `triangles.use`. ",
          "length(triangles) = ", length(triangles),
          ", length(Ys) = ", length(Ys), "."
        )
      )
    }

    if (length(triangles.wb) != length(Ys)) {
      stop(
        paste0(
          "disq_main(): aligned triangle list length does not match `Ys`. ",
          "length(triangles.wb) = ", length(triangles.wb),
          ", length(Ys) = ", length(Ys), "."
        )
      )
    }

    if (!is.null(fit_int$nonempty.index)) {
      cs.wb <- cs[fit_int$nonempty.index]
    } else if (length(cs) == length(Ys)) {
      cs.wb <- cs
    } else if (length(cs) == length(triangles.wb)) {
      cs.wb <- cs[seq_along(triangles.wb)]
    } else {
      stop(
        paste0(
          "disq_main(): cannot align `cs` with nonempty local subregions. ",
          "Please have `SQBiT_int()` return `nonempty.index`. ",
          "length(cs) = ", length(cs),
          ", length(Ys) = ", length(Ys),
          ", length(triangles.wb) = ", length(triangles.wb), "."
        )
      )
    }

    proj.H.wb <- fit_int$proj.H

    fit.vc.wb <- disq_vc_wb(
      tau = tau,
      lambdas = lambdas,
      theta.sm = thetas.sm,
      Ys = Ys,
      Cs = Cs,
      Xs = Xs,
      Ss = Ss,
      etas = eta.0,
      gammas = gammas,
      B.eva = B.eva,
      Br = Br,
      max.iter = max.iter,
      cs = cs.wb,
      eta.j = eta.j,
      d = d,
      r = r,
      level = level,
      Tr = Tr,
      V = V,
      triangles = triangles.wb,
      proj.H = proj.H.wb
    )

    beta.lb <- fit.vc.wb$beta.lb
    beta.ub <- fit.vc.wb$beta.ub
  }

  # ------------------------------------------------------------
  # Step 4. Select the master machine automatically if not provided
  # ------------------------------------------------------------
  if (is.null(master)) {
    master <- which.min(
      vapply(fit_int$resid, function(x) mean(rhotau(x, tau = tau)), numeric(1))
    )
  }

  # ------------------------------------------------------------
  # Step 5. Communication-based aggregation for constant coefficients
  # ------------------------------------------------------------
  fit_disq <- disq(
    tau = tau,
    m = m,
    maxTt = maxTt,
    c = c,
    Ys = Ys,
    Cs = Cs,
    Ws = Ws,
    gammas = gammas,
    Gs = Gs,
    eta = eta.0,
    alpha = alpha,
    beta1 = beta1,
    beta2 = beta2,
    master = master,
    level = level,
    method = method,
    stepMax = stepMax
  )

  end_time <- Sys.time()
  result_time <- end_time - start_time

  eta_init <- matrix(apply(eta.0, 2, mean), nrow = 1)
  eta_comm_trace <- fit_disq$trace[, -1, drop = FALSE]

  list(
    time = result_time,
    eta = fit_disq$eta[seq_len(q) + 1],
    eta.full = fit_disq$eta,
    betas = betas,
    thetas.sm = thetas.sm,
    beta.lb = beta.lb,
    beta.ub = beta.ub,
    fit.vc.wb = fit.vc.wb,
    cis = fit_disq$cis,
    ci.trace = fit_disq$ci.trace,
    eta.trace = rbind(
      eta_init[, -1, drop = FALSE],
      eta_comm_trace[-1, , drop = FALSE]
    ),
    fit.int = fit_int,
    fit.disq = fit_disq,
    master = master,
    m = m
  )
}
