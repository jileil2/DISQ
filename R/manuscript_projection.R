# Accumulate the weighted Gram in blocks to avoid a dense global design.
.disq_global_projection <- function(B, Q, X, C, w) {
  j <- ncol(Q)
  p <- ncol(X)
  n <- nrow(B)
  gram <- matrix(0, p * j, p * j)
  rhs <- matrix(0, p * j, ncol(C))
  columns <- split(seq_len(p * j), rep(seq_len(p), each = j))

  for (a in seq_len(p)) {
    aa <- columns[[a]]
    rhs[aa, ] <- as.matrix(Matrix::crossprod(
      Q, Matrix::crossprod(B, C * (w * X[, a]))
    )) / n
    for (b in seq_len(a)) {
      bb <- columns[[b]]
      raw <- Matrix::crossprod(B, Matrix::Diagonal(x = w * X[, a] * X[, b]) %*% B)
      block <- as.matrix(Matrix::crossprod(Q, raw %*% Q)) / n
      gram[aa, bb] <- block
      gram[bb, aa] <- t(block)
    }
  }

  scales <- sqrt(diag(gram))
  if (!all(is.finite(scales) & scales > 0)) stop("Nonpositive weighted Gram diagonal.")
  scaled_gram <- gram / outer(scales, scales)
  factor <- chol((scaled_gram + t(scaled_gram)) / 2)
  coefficient <- backsolve(factor, forwardsolve(t(factor), rhs / scales)) / scales
  prediction <- matrix(0, n, ncol(C))
  for (a in seq_len(p)) {
    spline <- Q %*% coefficient[columns[[a]], , drop = FALSE]
    prediction <- prediction + as.matrix(B %*% spline) * X[, a]
  }
  normal_error <- max(abs(gram %*% coefficient - rhs)) / max(1, max(abs(rhs)))
  if (!is.finite(normal_error) || normal_error >= 1e-8) {
    stop("Global weighted projection normal equations did not pass.")
  }

  list(D = C - prediction, P = t(coefficient), dimension = p * j,
       orthogonality_error = normal_error,
       normalized_gram_rcond = rcond(scaled_gram), condition = 1 / rcond(gram))
}
