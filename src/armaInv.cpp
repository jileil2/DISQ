// src/armaInv.cpp
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Matrix inverse (Armadillo)
//'
//' @param x A square numeric matrix.
//' @return The inverse of x.
//' @export
// [[Rcpp::export]]
arma::mat armaInv(const arma::mat& x) {
  // Basic safety check
  if (x.n_rows != x.n_cols) {
    Rcpp::stop("armaInv(): input matrix must be square.");
  }

  // Use solve instead of inv for better numerical behavior
  arma::mat I = arma::eye<arma::mat>(x.n_rows, x.n_cols);
  return arma::solve(x, I);
}