# ============================================================
# AInv(): stabilized matrix inverse via singular value decomposition
# ============================================================

#' Stabilized Matrix Inverse via Singular Value Decomposition
#'
#' \code{AInv()} computes a numerically stabilized inverse of a square
#' matrix using singular value decomposition (SVD).
#'
#' Given a matrix \eqn{A}, the function first computes its singular value
#' decomposition
#' \deqn{
#' A = U D V^\top,
#' }
#' where \eqn{D} is the diagonal matrix of singular values. The inverse is
#' then approximated by
#' \deqn{
#' A^{-1} \approx V D^{-1} U^\top.
#' }
#'
#' To avoid numerical failure when \eqn{A} is singular or nearly singular,
#' any zero singular values are replaced by a very small positive constant.
#' This yields a stabilized generalized inverse-like approximation that is
#' often useful in optimization and estimating-equation routines.
#'
#' @param A A numeric square matrix.
#' @param tol Small positive value used to replace zero singular values
#'   before inversion. The default is \code{1e-20}.
#'
#' @return A numeric matrix representing the stabilized inverse of \code{A}.
#'
#' @details
#' This function is intended for settings where a direct matrix inverse may
#' be unstable or undefined because of exact or near singularity. Rather
#' than failing when some singular values are zero, the function replaces
#' those values with a small constant before inversion.
#'
#' If
#' \deqn{
#' A = U \, \mathrm{diag}(d_1,\ldots,d_r) \, V^\top,
#' }
#' the function returns
#' \deqn{
#' V \, \mathrm{diag}(1/d_1,\ldots,1/d_r) \, U^\top,
#' }
#' after replacing any zero singular values by \code{tol}.
#'
#' This is not exactly the Moore-Penrose pseudoinverse, because zero
#' singular values are not dropped; instead they are regularized by a small
#' positive threshold.
#'
#' @examples
#' \donttest{
#' A <- matrix(c(2, 1, 1, 2), nrow = 2)
#' AInv(A)
#'
#' B <- matrix(c(1, 2, 2, 4), nrow = 2) # singular matrix
#' AInv(B)
#' }
#'
#' @export
AInv <- function(A, tol = 1e-20) {
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  }

  if (!is.numeric(A)) {
    stop("AInv(): 'A' must be a numeric matrix.")
  }

  if (nrow(A) != ncol(A)) {
    stop("AInv(): 'A' must be a square matrix.")
  }

  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("AInv(): 'tol' must be a positive finite scalar.")
  }

  svd.A <- svd(A)
  d <- svd.A$d
  u <- svd.A$u
  v <- svd.A$v

  d[d == 0] <- tol
  D.inv <- diag(1 / d, nrow = length(d), ncol = length(d))

  eigenMapMatMult(eigenMapMatMult(v, D.inv), t(u))
}
