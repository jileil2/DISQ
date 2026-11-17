# ============================================================
# kr(): row-wise Kronecker product
# ============================================================

#' Row-wise Kronecker Product
#'
#' Computes the row-wise Kronecker product (row-wise tensor product) of two matrices.
#' For each row \eqn{i}, this function returns \eqn{\mathrm{vec}(A_{i,\cdot} \otimes B_{i,\cdot})},
#' implemented as \code{as.vector(outer(A[i,], B[i,], "*"))}.
#'
#' This operator is used throughout the DISQ/SQBiT codebase to construct the design matrix
#' for varying-coefficient terms, e.g. \eqn{Z = \mathrm{kr}(X, BQ2)} where \eqn{BQ2 = B Q_2}.
#'
#' @param A Numeric matrix of dimension \eqn{n \times p_1}.
#' @param B Numeric matrix of dimension \eqn{n \times p_2}. Must have the same number of rows as \code{A}.
#' @param byrow Logical. Only \code{byrow = TRUE} is supported in DISQ.
#'
#' @return A numeric matrix of dimension \eqn{n \times (p_1 p_2)} whose \eqn{i}th row equals
#' the vectorized outer product of the \eqn{i}th rows of \code{A} and \code{B}.
#'
#' @examples
#' A <- matrix(1:6, nrow = 2, byrow = TRUE) # 2 x 3
#' B <- matrix(1:4, nrow = 2, byrow = TRUE) # 2 x 2
#' kr(A, B)
#'
#' # For i=1: outer(c(1,2,3), c(1,2)) -> c(1,2,2,4,3,6)
#'
#' @export
kr <- function(A, B, byrow = TRUE) {
  if (!isTRUE(byrow)) stop("kr(): only byrow=TRUE is supported in DISQ.")
  A <- as.matrix(A)
  B <- as.matrix(B)
  if (nrow(A) != nrow(B)) stop("kr(): A and B must have the same number of rows.")

  n <- nrow(A)
  p1 <- ncol(A)
  p2 <- ncol(B)

  out <- matrix(0, n, p1 * p2)
  for (i in seq_len(n)) out[i, ] <- as.vector(outer(A[i, ], B[i, ], "*"))
  out
}
