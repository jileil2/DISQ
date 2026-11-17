\name{SQBiT_int_metrics}
\alias{SQBiT_int_metrics}
\title{Local-fit residual diagnostics}
\description{An earlier diagnostic routine for distributed GACV, residual
spread, and Moran's I. It uses the original \code{tune.lambda} fitting helper.}
\usage{
SQBiT_int_metrics(tau, K1, maxTt = 100, Y, C, X, Q2,
  max.iter = 50, cs, eta.j = 1, d = 3, r = 1, Tr, V, S,
  triangles = FALSE, nlambda = 10, new.nlambda = 10,
  lambda.scale = 100, moran_k = 10, moran_style = "W",
  moran_mc = FALSE, moran_nsim = 199, center_resid = TRUE,
  sym_nb = FALSE)
}
\arguments{
\item{tau}{Quantile level.}
\item{K1}{Number of triangle groups for automatic partitioning.}
\item{maxTt, max.iter, eta.j, Q2}{Arguments retained by the earlier fitting interface.}
\item{Y}{Response vector.}
\item{C, X}{Constant and varying covariate matrices.}
\item{cs}{Local bandwidth multipliers.}
\item{d, r}{Spline degree and smoothness order.}
\item{Tr, V}{Triangle and vertex matrices.}
\item{S}{Observation coordinates.}
\item{triangles}{Triangle groups, or \code{TRUE} for automatic grouping.}
\item{nlambda, new.nlambda, lambda.scale}{Penalty search controls.}
\item{moran_k}{Number of nearest neighbors.}
\item{moran_style}{Spatial weights style passed to \code{spdep::nb2listw}.}
\item{moran_mc}{Whether to use a permutation test.}
\item{moran_nsim}{Number of permutations.}
\item{center_resid}{Whether to center residuals for Moran's I.}
\item{sym_nb}{Whether to symmetrize the neighbor graph.}
}
\value{A named numeric vector containing distributed GACV, residual spread, and Moran's I results.}
\keyword{internal}
