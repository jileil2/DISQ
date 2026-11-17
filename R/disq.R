# ============================================================
# disq(): Distributed aggregation for constant coefficient estimation
# ============================================================

#' Distributed Aggregation for Constant Coefficient Estimation
#'
#' Implements the multi-round communication algorithm for estimating
#' the constant coefficients in the decorrelated SQBiT (dSQBiT) model.
#' The procedure aggregates local estimators computed on multiple
#' machines and iteratively refines the global estimator using a
#' communication-efficient surrogate loss.
#'
#' At each communication round, local machines evaluate decorrelated
#' score functions at the current estimate and send the resulting
#' low-dimensional score vectors to the master machine. The master
#' aggregates these score vectors to approximate the global gradient
#' and updates the constant coefficient estimator by minimizing a
#' surrogate loss using either Barzilai-Borwein gradient descent or
#' the Adam optimizer.
#'
#' @param tau Quantile level.
#' @param m Number of local machines.
#' @param maxTt Maximum number of communication rounds.
#' @param c Positive constant used in the rule-of-thumb bandwidths for
#'   communication.
#' @param Ys List of response vectors for each machine.
#' @param Cs List of constant-covariate matrices for each machine.
#' @param Ws List of basis matrices associated with varying coefficients.
#' @param gammas List of estimated varying-coefficient parameters from
#'   the initial SQBiT step.
#' @param Gs List of estimated projection matrices for each machine.
#' @param eta Matrix of local initial estimators, including the intercept.
#' @param alpha Learning rate for the Adam optimizer.
#' @param beta1 Exponential decay rate for the first moment in Adam.
#' @param beta2 Exponential decay rate for the second moment in Adam.
#' @param master Index of the master machine.
#' @param level Significance level for confidence intervals.
#' @param method Optimization method used in each communication round.
#'   Supported options are `"GD-BB"` and `"Adam"`.
#' @param stepMax Maximum step size passed to the inner optimizer.
#'
#' @return A list with components:
#' \describe{
#'   \item{eta}{Final aggregated estimator of the constant coefficients,
#'   including the intercept.}
#'   \item{trace.sc}{Trace of aggregated gradients across communication rounds,
#'   excluding the intercept.}
#'   \item{time}{Total running time.}
#'   \item{trace}{Trace of coefficient estimates across communication rounds.}
#'   \item{ci.trace}{Stored confidence intervals at the final estimator.}
#'   \item{resid}{Residuals at the master machine evaluated at the final estimator.}
#'   \item{nCt}{Number of communication rounds completed.}
#'   \item{gH}{Estimated Hessian-type matrix used for inference.}
#'   \item{J1.inv}{Inverse of the estimated information matrix.}
#'   \item{J1}{Estimated information matrix.}
#'   \item{cis}{Confidence intervals for the constant coefficients,
#'   excluding the intercept.}
#'   \item{asd}{Estimated asymptotic standard deviations.}
#'   \item{trace.ns}{Mean absolute aggregated gradient by round,
#'   excluding the intercept.}
#'   \item{ns}{Vector of local sample sizes computed from `Ys`.}
#' }
#'
#' @export
disq <- function(tau, m, maxTt = 100, c = 2.5,
                 Ys, Cs, Ws, gammas, Gs, eta,
                 alpha = 1e-5, beta1 = 0.9, beta2 = 0.999,
                 master, level = 0.05, method = "GD-BB",
                 stepMax = 1e-3) {
  start_time <- Sys.time()

  # ---------------------------
  # basic coercions
  # ---------------------------
  Ys <- lapply(Ys, as.numeric)
  Cs <- lapply(Cs, function(x) as.matrix(x))
  Ws <- lapply(Ws, function(x) as.matrix(x))
  Gs <- lapply(Gs, function(x) as.matrix(x))
  gammas <- lapply(gammas, function(x) matrix(as.numeric(x), ncol = 1))
  eta <- as.matrix(eta)

  # local sample sizes
  ns <- vapply(Ys, length, integer(1))

  # basic info
  N <- sum(ns)
  ws <- ns / N
  n.m <- ns[master]

  # number of constant coefficients excluding intercept
  q <- ncol(eta) - 1L
  if (q <= 0L) {
    stop("disq(): eta must include an intercept and at least one constant coefficient.")
  }

  # initial estimates (includes intercept)
  eta.t <- matrix(apply(eta, 2, mean), ncol = 1)

  # communication bandwidths
  h <- c * ((q + log(N)) / N)^(1 / 3)
  b <- c * ((q + log(n.m)) / n.m)^(1 / 3)

  score <- function(G, Y, gamma, eta, W, C, h, tau) {
    G <- as.matrix(G)
    W <- as.matrix(W)
    C <- as.matrix(C)
    gamma <- matrix(as.numeric(gamma), ncol = 1)
    eta <- matrix(as.numeric(eta), ncol = 1)

    n <- nrow(W)
    Y.star <- as.numeric(Y) - as.numeric(eigenMapMatMult(W, gamma))
    r <- Y.star - as.numeric(eigenMapMatMult(C, eta))
    C.star <- as.matrix(eigenMapMatMult(G, C))
    sc <- matrix(punif(-r / h, -1, 1) - tau, ncol = 1)

    # returns a column vector of length q+1
    out <- as.matrix(eigenMapMatMult(t(C.star), sc)) / n
    matrix(as.numeric(out), ncol = 1)
  }

  # master objects
  C1 <- Cs[[master]]
  Y1 <- Ys[[master]]
  gamma1 <- gammas[[master]]
  W1 <- Ws[[master]]

  # projected constant covariates
  Cs.proj <- vector("list", m)
  for (k in seq_len(m)) {
    Cs.proj[[k]] <- as.matrix(eigenMapMatMult(Gs[[k]], Cs[[k]]))
  }
  C1.proj <- Cs.proj[[master]]

  # traces
  trace <- matrix(as.numeric(eta.t), nrow = 1)
  gt <- numeric(0)
  trace.sc <- matrix(0, nrow = 1, ncol = q)
  ci.trace <- vector("list", 0)

  for (t in seq_len(maxTt)) {
    grad.1.b <- score(
      Y = Ys[[master]],
      C = cbind(1, Cs[[master]]),
      W = Ws[[master]],
      G = Gs[[master]],
      tau = tau,
      eta = eta.t,
      gamma = gammas[[master]],
      h = b
    )

    grads.t <- ws[master] * score(
      Y = Ys[[master]],
      C = cbind(1, Cs[[master]]),
      W = Ws[[master]],
      G = Gs[[master]],
      tau = tau,
      eta = eta.t,
      gamma = gammas[[master]],
      h = h
    )

    for (a in seq_len(m)[-master]) {
      grads.t <- grads.t + ws[a] * score(
        Y = Ys[[a]],
        C = cbind(1, Cs[[a]]),
        W = Ws[[a]],
        G = Gs[[a]],
        tau = tau,
        eta = eta.t,
        gamma = gammas[[a]],
        h = h
      )
    }

    grads.t <- matrix(as.numeric(grads.t), ncol = 1)
    grad.cc <- as.numeric(grads.t[-1, 1])

    if (length(grad.cc) != q) {
      stop(
        sprintf(
          "disq(): gradient dimension mismatch. length(grads.t[-1]) = %d but expected %d (= ncol(eta) - 1).",
          length(grad.cc), q
        )
      )
    }

    gt <- c(gt, mean(abs(grad.cc)))
    trace.sc <- rbind(trace.sc, matrix(grad.cc, nrow = 1, ncol = q))

    stop_now <- (t > 1) && (gt[t] < min(1e-4, 1 / sqrt(N)))
    diverge <- (t > 1) && (gt[t] > gt[t - 1])

    if (diverge || stop_now) {
      break
    }

    grad.diff <- grad.1.b - grads.t

    if (method == "GD-BB") {
      eta.t <- smqr.eta.dc(
        h = b,
        Y = Ys[[master]],
        C = cbind(1, Cs[[master]]),
        W = Ws[[master]],
        gamma = gammas[[master]],
        G = Gs[[master]],
        initial = eta.t,
        tau = tau,
        grad = grad.diff,
        epsilon = 1e-4,
        stepMax = stepMax
      )
      eta.t <- matrix(as.numeric(eta.t), ncol = 1)
    } else if (method == "Adam") {
      eta.t <- smqr.eta.dc_adam(
        h = b,
        Y = Ys[[master]],
        C = cbind(1, Cs[[master]]),
        W = Ws[[master]],
        gamma = gammas[[master]],
        G = Gs[[master]],
        alpha = alpha,
        beta1 = beta1,
        beta2 = beta2,
        initial = eta.t,
        tau = tau,
        grad = grad.diff,
        epsilon = 1e-4,
        stepMax = stepMax
      )
      eta.t <- matrix(as.numeric(eta.t), ncol = 1)
    } else {
      stop("disq(): method must be either 'GD-BB' or 'Adam'.")
    }

    trace <- rbind(trace, matrix(as.numeric(eta.t), nrow = 1))
  }

  # inference computed once at the final estimator
  resid <- as.numeric(
    Y1 - as.numeric(eigenMapMatMult(W1, gamma1)) -
      as.numeric(eigenMapMatMult(cbind(1, C1), eta.t))
  )

  bounds <- quantile(resid, c(0.025, 0.975))
  rr <- bounds[2] - bounds[1]
  b.ve <- 0.45 * rr^2 * (q * log(n.m) / n.m)^(1 / 3)
  kd <- dnorm(resid / b.ve) / b.ve
  Sigma.h <- diag(kd)

  C1.w.prod <- as.matrix(eigenMapMatMult(Sigma.h, C1))
  C1.w.prod <- as.matrix(eigenMapMatMult(t(C1.proj), C1.w.prod))
  gH <- C1.w.prod / n.m

  J1 <- gH
  J1.inv <- AInv(J1)

  Xi.N.b <- as.matrix(eigenMapMatMult(J1.inv, gH))
  Xi.N.b <- as.matrix(eigenMapMatMult(Xi.N.b, J1.inv))
  asd <- sqrt(diag(Xi.N.b))
  Xi.N.b <- tau * (1 - tau) * Xi.N.b

  cis <- matrix(NA_real_, nrow = ncol(C1), ncol = 2)
  for (k in seq_len(ncol(C1))) {
    cis[k, ] <- eta.t[k + 1, 1] +
      c(-1, 1) * qnorm(1 - level / 2) * sqrt(Xi.N.b[k, k] / N)
  }
  ci.trace[[1]] <- cis

  end_time <- Sys.time()
  result_time <- end_time - start_time

  list(
    eta = eta.t,
    trace.sc = trace.sc,
    time = result_time,
    trace = trace,
    ci.trace = ci.trace,
    resid = resid,
    nCt = nrow(trace),
    gH = gH,
    J1.inv = J1.inv,
    J1 = J1,
    cis = cis,
    asd = asd,
    trace.ns = gt,
    ns = ns
  )
}
