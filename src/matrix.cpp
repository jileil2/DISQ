// Matrix operations used by local fitting and communication.
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Kronecker product (Armadillo)
//'
//' @param a numeric matrix
//' @param b numeric matrix
//' @return kron(a, b)
//' @export
// [[Rcpp::export]]
arma::mat cppkp(const arma::mat& a, const arma::mat& b) {
  return arma::kron(a, b);
}

//' Matrix multiplication (Armadillo)
//'
//' @param A numeric matrix
//' @param B numeric matrix
//' @return A %*% B
//' @export
// [[Rcpp::export]]
arma::mat eigenMatMult(const arma::mat& A, const arma::mat& B) {
  return A * B;
}

//' Matrix multiplication (alias; Armadillo)
//'
//' Retains the original matrix-multiplication interface.
//'
//' @param A numeric matrix
//' @param B numeric matrix
//' @return A %*% B
//' @export
// [[Rcpp::export]]
arma::mat eigenMapMatMult(const arma::mat& A, const arma::mat& B) {
  return A * B;
}

//' Matrix multiplication with n_cores (alias; Armadillo)
//'
//' Note: Armadillo/BLAS threading is typically controlled via environment
//' variables; this function keeps the same interface but does not force
//' the thread count.
//'
//' @param A numeric matrix
//' @param B numeric matrix
//' @param n_cores number of cores (unused here; kept for compatibility)
//' @return A %*% B
//' @export
// [[Rcpp::export]]
arma::mat eigenMapMatMult2(const arma::mat& A, const arma::mat& B, const int n_cores) {
  (void)n_cores; // suppress unused warning
  return A * B;
}

//' Indicator for negative
//'
//' @param x number
//' @return 1 if x < 0 else 0
//' @export
// [[Rcpp::export]]
int neg(const double x) {
  return (x < 0.0) ? 1 : 0;
}

//' Sum of a numeric vector
//'
//' @param x numeric vector
//' @return sum(x)
//' @export
// [[Rcpp::export]]
double sum_cpp(const Rcpp::NumericVector& x) {
  double s = 0.0;
  const R_xlen_t n = x.size();
  for (R_xlen_t i = 0; i < n; ++i) s += x[i];
  return s;
}

//' Sum of check loss rho_tau over a vector
//'
//' rho_tau(u) = u * (tau - I(u < 0))
//'
//' @param x numeric vector (residuals)
//' @param tau quantile level in (0,1)
//' @return sum_i rho_tau(x_i)
//' @export
// [[Rcpp::export]]
double sum_rhotau(const Rcpp::NumericVector& x, const double tau) {
  if (!(tau > 0.0 && tau < 1.0)) {
    Rcpp::stop("sum_rhotau(): tau must be in (0, 1).");
  }

  double s = 0.0;
  const R_xlen_t n = x.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    const double u = x[i];
    s += u * (tau - (u < 0.0 ? 1.0 : 0.0));
  }
  return s;
}

//' Pairwise Euclidean distance between rows of X (Armadillo)
//'
//' @param X numeric matrix, rows are points
//' @return n x n distance matrix
//' @export
// [[Rcpp::export]]
arma::mat arma_dist(const arma::mat& X) {
  const arma::uword n = X.n_rows;
  arma::mat D(n, n, arma::fill::zeros);

  for (arma::uword i = 0; i < n; ++i) {
    D(i, i) = 0.0;
    for (arma::uword j = i + 1; j < n; ++j) {
      const double d = arma::norm(X.row(i) - X.row(j), 2);
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}