# Implicit master updates with a fixed global decorrelating projection.
.disq_communicate <- function(Y, C, S, nuisance, D, index, ns, master, h, b,
                              tau, max_rounds = 15L, tol = 1e-7,
                              vcov_radius, eta_initial, level = 0.05,
                              kernel = "norm") {
  start <- proc.time()[["elapsed"]]
  Y <- as.numeric(Y)
  C <- as.matrix(C)
  S <- as.matrix(S)
  nuisance <- as.numeric(nuisance)
  D <- as.matrix(D)
  eta_initial <- as.numeric(eta_initial)
  n <- length(Y)
  q <- ncol(C)
  k <- length(index)
  ids <- as.integer(unlist(index, use.names = FALSE))
  if (n < 1L || q < 1L || nrow(C) != n || nrow(S) != n || ncol(S) != 2L ||
    length(nuisance) != n || !identical(dim(D), dim(C)) ||
    length(eta_initial) != q || any(!is.finite(c(Y, C, S, nuisance, D, eta_initial))) ||
    k < 1L || length(ns) != k || !identical(as.integer(ns), as.integer(lengths(index))) ||
    length(ids) != n || anyDuplicated(ids) || !identical(sort(ids), seq_len(n))) {
    stop("Communication inputs must be finite, use original observation order, and form an exact partition.")
  }
  if (length(master) != 1L || !is.finite(master) || master != floor(master) ||
    master < 1L || master > k || any(ns <= 0)) stop("Invalid master or empty worker.")
  if (length(max_rounds) != 1L || !is.finite(max_rounds) ||
    max_rounds < 0 || max_rounds != floor(max_rounds) ||
    length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
    length(h) != 1L || !is.finite(h) || h <= 0 ||
    length(b) != 1L || !is.finite(b) || b <= 0 ||
    length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1 ||
    length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("Invalid communication bandwidths, tolerances, round count or quantile level.")
  }
  kernel <- match.arg(kernel, c("norm", "unif"))
  pair_kernel <- .disq_pair_kernel(S, vcov_radius)
  worker_score <- function(worker, eta, bandwidth) {
    i <- index[[worker]]
    .disq_score(eta, Y[i], C[i, , drop = FALSE], nuisance[i],
      D[i, , drop = FALSE], bandwidth, tau,
      kernel = kernel
    )
  }
  global_score <- function(eta) {
    Reduce(`+`, lapply(seq_len(k), function(worker)
      (ns[worker] / n) * worker_score(worker, eta, h)))
  }
  eta <- eta_initial
  trace <- matrix(eta, nrow = 1L)
  score_trace <- matrix(global_score(eta), nrow = 1L)
  root_diagnostics <- vector("list", max_rounds)
  failed_round <- NA_integer_
  failure <- NULL
  for (round in seq_len(max_rounds)) {
    current_score <- as.numeric(score_trace[nrow(score_trace), ])
    offset <- worker_score(master, eta, b) - current_score
    i <- index[[master]]
    root <- tryCatch(.disq_solve_score(
      eta, Y[i], C[i, , drop = FALSE], nuisance[i], D[i, , drop = FALSE],
      b, tau,
      offset = offset, max.iter = 200L,
      tolerance = min(tol / 10, 1e-9), kernel = kernel
    ), error = function(e) e)
    if (inherits(root, "error")) {
      failure <- conditionMessage(root)
      root_diagnostics[[round]] <- list(round = round, converged = FALSE, error = failure)
      failed_round <- as.integer(round)
      break
    }
    root_diagnostics[[round]] <- c(
      list(round = round),
      root[c("converged", "score_norm", "iterations", "method", "status")]
    )
    if (!isTRUE(root$converged) || any(!is.finite(root$eta))) {
      failure <- paste("Master root did not converge:", root$status)
      failed_round <- as.integer(round)
      break
    }
    eta <- as.numeric(root$eta)
    new_score <- global_score(eta)
    if (any(!is.finite(new_score))) {
      failure <- "Nonfinite global score after master root."
      failed_round <- as.integer(round)
      break
    }
    # Every stored row is an actual evaluation of the implicit update. There
    # is no early stopping at tol and no carrying endpoints into missing rows.
    trace <- rbind(trace, eta)
    score_trace <- rbind(score_trace, new_score)
  }
  rounds <- seq_len(nrow(trace)) - 1L
  root_diagnostics <- Filter(Negate(is.null), root_diagnostics)
  covariance <- lapply(seq_len(nrow(trace)), function(i)
    .disq_covariance(as.numeric(trace[i, ]), Y, C, nuisance, D, h, tau,
      pair_kernel,
      level = level, kernel = kernel
    ))
  final <- covariance[[length(covariance)]]
  score_norm <- apply(abs(score_trace), 1L, max)
  mean_absolute_score <- rowMeans(abs(score_trace))
  list(
    eta = as.numeric(tail(trace, 1L)), trace = trace, eta.trace = trace,
    rounds = as.integer(rounds), score_trace = score_trace, trace.sc = score_trace,
    score_norm = tail(score_norm, 1L), score_norm_trace = score_norm,
    mean_absolute_score = mean_absolute_score,
    manuscript_score_met = mean_absolute_score < 1e-4,
    tolerance_reached = score_norm <= tol,
    completed = is.null(failure) && length(rounds) == max_rounds + 1L,
    converged = is.null(failure) && tail(score_norm, 1L) <= tol,
    failed_round = failed_round, error = failure,
    missing_rounds = setdiff(seq.int(0L, as.integer(max_rounds)), rounds),
    root_diagnostics = root_diagnostics,
    covariance = covariance,
    V = final$V, J = final$J, Gamma = final$Gamma, SE = final$se, CI = final$ci,
    V_trace = lapply(covariance, `[[`, "V"),
    J_trace = lapply(covariance, `[[`, "J"),
    Gamma_trace = lapply(covariance, `[[`, "Gamma"),
    SE_trace = do.call(rbind, lapply(covariance, `[[`, "se")),
    CI_trace = lapply(covariance, `[[`, "ci"),
    covariance_valid = vapply(covariance, `[[`, logical(1), "valid"),
    h = h, b = b, tau = tau, master = master, ns = ns,
    vcov_radius = vcov_radius, kernel = kernel,
    runtime_seconds = proc.time()[["elapsed"]] - start
  )
}
