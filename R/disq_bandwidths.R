disq_bandwidths <- function(N, J, n_master, J_master, c = 0.20) {
  counts <- list(N = N, J = J, n_master = n_master, J_master = J_master)
  valid <- vapply(counts, function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  }, logical(1))
  if (!all(valid)) {
    stop("Supply a positive integer for ", paste(names(counts)[!valid], collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (n_master > N) stop("n_master cannot exceed N.", call. = FALSE)
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c) || c <= 0) {
    stop("c must be a positive finite number.", call. = FALSE)
  }

  log_size <- log(2 + N)^2
  base::c(
    h = c * (J * log_size / N)^(1 / 3),
    b = c * (J_master * log_size / n_master)^(1 / 3)
  )
}
