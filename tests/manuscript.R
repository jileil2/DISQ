library(DISQ)

stopifnot(isTRUE(all.equal(
  disq_bandwidths(82469, 626, 12490, 86),
  c(h = 0.19816516573715062, b = 0.19182775795481893),
  tolerance = 1e-14
)))

set.seed(1328)
n <- 800L
index <- list(1:400, 401:800)
B <- Matrix::sparseMatrix(i = seq_len(n), j = rep(1:2, each = 400), x = 1)
Q2 <- matrix(1 / sqrt(2), 2, 1)
X <- matrix(1, n, 1)
C <- matrix(rnorm(n), n, 1)
S <- cbind(seq_len(n), 0)
Y <- 0.2 + 0.7 * C[, 1L] + rnorm(n, sd = 0.7)
local <- list(
  list(columns = 1, Q2 = matrix(1), theta = matrix(0.18), eta = 0.68, phi = 0.35),
  list(columns = 2, Q2 = matrix(1), theta = matrix(0.23), eta = 0.72, phi = 0.35)
)
fit <- disq_from_local_fits(Y, C, X, S, B, Q2, index, local, tau = 0.5, radius = 0.25)
stopifnot(fit$completed, fit$converged, identical(fit$communication$rounds, 0:15),
          fit$covariance$valid, length(fit$eta) == 1,
          identical(fit$communication$mean_absolute_score,
                    rowMeans(abs(fit$communication$score_trace))),
          isTRUE(tail(fit$communication$manuscript_score_met, 1L)),
          max(abs(fit$initialization$eta_mean - colMeans(fit$initialization$eta_star))) < 1e-14,
          max(abs(fit$beta - 0.205)) < 1e-14)

h <- fit$bandwidths[["h"]]
nuisance <- rowSums(X * fit$beta)
D <- fit$projection$D
weight <- dnorm((Y - nuisance - C[, 1] * fit$initialization$eta_mean) / h) / h
stopifnot(max(abs(crossprod(as.matrix(B %*% Q2), D * weight))) < 1e-10)
final_residual <- Y - nuisance - C[, 1] * fit$eta
score <- D * (pnorm(-final_residual / h) - .5)
jacobian <- crossprod(D, C * dnorm(final_residual / h) / h) / n
sandwich <- solve(jacobian) %*% crossprod(score) %*% t(solve(jacobian)) / n^2
stopifnot(max(abs(colMeans(score))) < 1e-7,
          isTRUE(all.equal(fit$covariance$V, sandwich, tolerance = 1e-12)))

all_varying <- disq_from_local_fits(Y, matrix(numeric(0), n, 0), X, S, B, Q2,
                                   index, local, tau = 0.5, radius = 0.25)
stopifnot(all_varying$completed, length(all_varying$eta) == 0L,
          identical(all_varying$fitted, rowSums(X * fit$beta)))
