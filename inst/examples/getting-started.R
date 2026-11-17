library(DISQ)

# Simulate a response on a square split into two triangular subregions.
set.seed(1328)
n <- 1600L
S <- cbind(runif(n), runif(n))
V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
Tr <- rbind(c(1, 2, 3), c(1, 3, 4))
index <- list(which(S[, 2] < S[, 1]), which(S[, 2] >= S[, 1]))
C <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "constant"))
X <- matrix(1, n, 1)
Y <- 0.2 + 0.3 * S[, 1] - 0.15 * S[, 2] + 0.7 * C[, 1] + rnorm(n, sd = 0.7)
tau <- 0.5

# Fit the varying intercept and constant slope separately on each worker.
global_basis <- BPST::basis(V, Tr, d = 3, r = 1, Z = S)
B <- global_basis$B
Q2 <- as.matrix(global_basis$Q2)
columns_per_triangle <- choose(3 + 2, 2)

local <- lapply(seq_along(index), function(k) {
  rows <- index[[k]]
  basis <- BPST::basis(V, Tr[k, , drop = FALSE], d = 3, r = 1, Z = S[rows, ])
  # A single triangle has no cross-triangle smoothness constraints.
  Q2_local <- diag(ncol(basis$B))
  penalty <- crossprod(Q2_local, basis$K %*% Q2_local)
  phi <- 7 * tau * (1 - tau) *
    ((ncol(C) + ncol(Q2_local) * ncol(X) + log(length(rows))) / length(rows))^(2 / 5)

  fit <- SQBiT::SQBiT_gacv(
    tau = tau, Y = Y[rows], C = C[rows, , drop = FALSE],
    X = X[rows, , drop = FALSE], B = as.matrix(basis$B),
    Q2 = Q2_local, P = as.matrix(penalty), h = phi, max.iter = 2000,
    eta.j1 = 0.32, eta.j2 = 0.32
  )
  list(
    columns = (k - 1L) * columns_per_triangle + seq_len(columns_per_triangle),
    Q2 = Q2_local, theta = as.matrix(fit$theta), eta = fit$eta, phi = phi
  )
})

# Combine the local fits and retain the complete communication trace.
fit <- disq_from_local_fits(
  Y, C, X, S, B, Q2, index, local,
  tau = tau,
  radius = 0.035, c = 0.20, rounds = 15, level = 0.01
)
fit$eta
fit$covariance$ci
fit$communication$mean_absolute_score

# The response in this example is untransformed, so fitted quantiles remain
# on its original scale. For a logit response, use plogis(fit$fitted).
