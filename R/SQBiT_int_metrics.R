# ============================================================
# SQBiT_int_metrics(): return GACV, GSD, Moran's I, and p-value
# (kNN weights built from spatial coordinates S)
# ============================================================

SQBiT_int_metrics <- function(tau, K1, maxTt = 100,
                              Y, C, X, Q2, max.iter = 50,
                              cs, eta.j = 1, d = 3, r = 1,
                              Tr, V, S, triangles = FALSE,
                              nlambda = 10, new.nlambda = 10, lambda.scale = 100,
                              # ---- Moran's I options ----
                              moran_k = 10,
                              moran_style = "W", # "W" row-standardized
                              moran_mc = FALSE, # TRUE -> moran.mc permutation test
                              moran_nsim = 199,
                              center_resid = TRUE, # center residuals before Moran
                              sym_nb = FALSE # TRUE -> symmetrize kNN graph
) {
  # ---- dependencies ----
  if (!requireNamespace("spdep", quietly = TRUE)) stop("Please install 'spdep'.")
  if (!requireNamespace("RANN", quietly = TRUE)) stop("Please install 'RANN'.")

  # ---- basic info ----
  N <- length(Y)

  # ---- split into triangles (if triangles is TRUE) ----
  if (is.logical(triangles) && triangles == TRUE) {
    triangles <- list()
    nTr <- nrow(Tr)
    indices <- 1:nTr
    nT <- ceiling(nTr / K1)

    for (k in 1:K1) {
      v.1 <- apply(V[Tr[indices[1], ], , drop = FALSE], 2, mean)

      ni <- length(indices) - 1
      dist <- numeric(ni)

      for (i in 1:ni) {
        v.i <- apply(V[Tr[indices[i + 1], ], , drop = FALSE], 2, mean)
        dist[i] <- norm(matrix(v.1) - matrix(v.i), "f")
      }

      take <- min(nT - 1, length(dist))
      cand.k <- c(indices[1], indices[order(dist)[1:take] + 1])

      triangles[[k]] <- as.numeric(cand.k)
      indices <- indices[which(!(indices %in% cand.k))]
    }

    triangles[[K1]] <- triangles[[K1]][1:(nTr - (K1 - 1) * nT)]
  }

  if (!is.list(triangles)) stop("'triangles' must be TRUE or a list of triangle-index sets.")

  # ---- keep only non-empty triangle groups ----
  nt <- length(triangles)
  indic <- integer(nt)

  for (a in 1:nt) {
    ind.a <- inVT(
      V0 = V, Tr0 = Tr[triangles[[a]], , drop = FALSE],
      xx = S[, 1], yy = S[, 2]
    )
    indic[a] <- as.integer(length(ind.a$ind.inside) > 0)
  }
  triangles <- triangles[which(indic == 1)]
  nt <- length(triangles)

  # Worker sizes and master index.
  ns <- numeric(nt)
  for (a in 1:nt) {
    ind.a <- inVT(
      V0 = V, Tr0 = Tr[triangles[[a]], , drop = FALSE],
      xx = S[, 1], yy = S[, 2]
    )
    ns[a] <- length(ind.a$ind.inside)
  }
  master <- which.max(ns)
  nq <- (d + 1) * (d + 2) / 2

  # ---- containers ----
  resid_list <- vector("list", nt)
  dfvs <- numeric(nt)
  mads <- numeric(nt)

  # For global residual alignment:
  resid_all <- rep(NA_real_, N)
  resid_cnt <- integer(N)

  # ---- loop over subregions / triangle groups ----
  for (a in 1:nt) {
    ind.a <- inVT(
      V0 = V, Tr0 = Tr[triangles[[a]], , drop = FALSE],
      xx = S[, 1], yy = S[, 2]
    )$ind.inside

    if (length(ind.a) == 0) next

    Y.a <- Y[ind.a]
    X.a <- X[ind.a, , drop = FALSE]
    S.a <- S[ind.a, , drop = FALSE]
    C.a <- as.matrix(C[ind.a, , drop = FALSE])
    n.a <- length(Y.a)

    # basis and penalty
    B0.sub <- basis(
      V = V, Tr = Tr[triangles[[a]], , drop = FALSE],
      d = d, r = r, Z = S.a
    )
    Q2.sub <- B0.sub$Q2
    B.sub <- B0.sub$B
    BQ2.sub <- as.matrix(B.sub %*% Q2.sub)

    K.sub <- B0.sub$K
    P.sub <- t(Q2.sub) %*% K.sub %*% Q2.sub

    # design for varying coefficients / spatial basis
    W.a <- as.matrix(kr(as.matrix(X.a), BQ2.sub, byrow = TRUE))

    # local bandwidth
    h <- cs[a] * tau * (1 - tau) *
      ((ncol(C.a) + ncol(Q2.sub) * ncol(X) + log(n.a)) / n.a)^(2 / 5)

    # tune smoothing parameter lambda by GACV (your routine)
    fit.a <- tune.lambda(
      Y = Y.a, X = as.matrix(X.a), C = C.a, max.iter = 100,
      P = P.sub, B = B.sub, Q2 = Q2.sub, var.j.2 = FALSE,
      kernel = "unif", h = h, tau = tau, lambda_start = 1e-4,
      cri = "gacv", eta.j1 = .32, eta.j2 = .32,
      nlambda = nlambda, new.nlambda = new.nlambda, lambda.scale = lambda.scale
    )

    dfvs[a] <- fit.a$df
    eta.0 <- fit.a$eta

    # residuals in this subregion
    resid.a <- as.numeric(
      Y.a -
        eigenMapMatMult(W.a, fit.a$gamma[-(1:ncol(C))]) -
        eigenMapMatMult(C.a, matrix(eta.0))
    )

    resid_list[[a]] <- resid.a
    mads[a] <- mad(resid.a)

    # ---- align residuals to full index space ----
    ok <- !is.na(resid.a)
    idxi <- ind.a[ok]

    resid_all[idxi] <- ifelse(is.na(resid_all[idxi]), 0, resid_all[idxi]) + resid.a[ok]
    resid_cnt[idxi] <- resid_cnt[idxi] + 1L
  }

  # ---- finalize global residual vector (average if overlaps) ----
  has <- resid_cnt > 0L
  resid_all[has] <- resid_all[has] / resid_cnt[has]

  if (any(!has)) {
    warning(sprintf(
      "There are %d points with no residual assigned; Moran's I will omit them.",
      sum(!has)
    ))
  }

  # ---- compute GACV ----
  resid_vec <- as.vector(unlist(resid_list))
  loss <- rhotau(u = resid_vec, tau = tau)
  gacv <- sum(loss) / (N - sum(dfvs))

  # ---- compute GSD ----
  w <- dfvs / sum(dfvs)
  mu <- sum(w * log(mads))
  H <- sum(w * (log(mads) - mu)^2)
  GSD <- exp(sqrt(H))

  # ---- Moran's I on residuals ----
  resid_use <- resid_all[has]
  S_use <- S[has, , drop = FALSE]
  n_use <- length(resid_use)

  if (center_resid) resid_use <- resid_use - mean(resid_use)

  # Build kNN neighbors via RANN (fast for big n)
  k <- moran_k
  nn <- RANN::nn2(S_use, k = k + 1) # includes self
  idx <- nn$nn.idx[, -1, drop = FALSE] # remove self

  # Build nb object from neighbor indices
  nb <- lapply(seq_len(n_use), function(i) idx[i, ])
  class(nb) <- "nb"
  attr(nb, "region.id") <- as.character(seq_len(n_use))
  attr(nb, "type") <- "knn"
  attr(nb, "sym") <- FALSE

  if (sym_nb) nb <- spdep::make.sym.nb(nb)

  lw <- spdep::nb2listw(nb, style = moran_style, zero.policy = TRUE)

  # Moran's I + p-value
  if (!moran_mc) {
    mt <- spdep::moran.test(resid_use, lw, zero.policy = TRUE)
    Moran_I <- unname(mt$estimate["Moran I statistic"])
    Moran_p <- unname(mt$p.value)
  } else {
    mc <- spdep::moran.mc(resid_use, lw, nsim = moran_nsim, zero.policy = TRUE)
    Moran_I <- unname(mc$statistic)
    # moran.mc returns p.value (Monte Carlo)
    Moran_p <- unname(mc$p.value)
  }

  return(c(gacv = gacv, GSD = GSD, Moran_I = Moran_I, Moran_p = Moran_p))
}
