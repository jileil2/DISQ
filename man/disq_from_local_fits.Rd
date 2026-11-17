\name{disq_from_local_fits}
\alias{disq_from_local_fits}
\title{Combine local SQBiT fits by fixed-projection communication}
\usage{
disq_from_local_fits(Y, C, X, S, B, Q2, index, local, tau,
  radius, master = which.max(lengths(index)), c = 0.20,
  rounds = 15L, level = 0.05, score_tolerance = 1e-7,
  B_eval = NULL)
}
\arguments{
  \item{Y}{Response in original observation order.}
  \item{C, X}{Matrices for the constant and varying terms. Include a column
  of ones in \code{X} for a varying intercept. Use a zero-column matrix for
  \code{C} when all coefficients vary.}
  \item{S}{Two-column spatial coordinates in the units of the covariance cutoff.}
  \item{B, Q2}{Global raw Bernstein design and orthonormal basis for its
  smoothness constraints. Their product is the constrained spline design.}
  \item{index}{List of observation row indices for the workers.}
  \item{local}{One list per worker, in the order of \code{index}. Each entry
  contains \code{columns} indexing its raw columns in \code{B}, its
  orthonormal constraint basis \code{Q2}, the fitted raw spline coefficient
  matrix \code{theta}, the local constant estimates \code{eta}, and the
  positive Gaussian bandwidth \code{phi}. Worker raw columns partition
  the global raw columns. The rows of \code{theta} follow \code{columns}.}
  \item{tau}{Quantile level, strictly between zero and one.}
  \item{radius}{Positive cutoff for the finite-range spatial covariance.}
  \item{master}{Worker used for the implicit update. Defaults to the largest worker.}
  \item{c}{Common communication-bandwidth multiplier.}
  \item{rounds}{Number of communication updates requested.}
  \item{level}{Significance level for constant-coefficient intervals.}
  \item{score_tolerance}{Global maximum absolute score used to report convergence.}
  \item{B_eval}{Optional raw global spline basis at prediction sites.}
}
\description{Performs the local decorrelated initialization, smoothness correction,
global density-weighted projection, fixed-projection communication, and
finite-range spatial covariance calculation.}
\details{Supply the local SQBiT fits at the selected allocation and quantile.
The common projection uses a Gaussian kernel and remains fixed during
communication. The function records each actual update and retains a partial
trace when a master root fails. Only the requested endpoint supplies the
returned constant estimates and their covariance. No wild-bootstrap calculations
are performed by this function.}
\value{A list with the final constant estimates \code{eta}, constrained
raw spline coefficients \code{theta}, varying estimates \code{beta}, optional
evaluations \code{beta_eval}, and fitted conditional quantiles \code{fitted}.
For models with constants it also contains local root diagnostics, the projection,
bandwidths, communication trace and endpoint spatial covariance. The
\code{completed} flag records availability of the requested endpoint.
The \code{converged} flag reports the endpoint maximum absolute score check.
The communication trace also reports the mean absolute score and whether
it is below the manuscript threshold of \eqn{10^{-4}} at each iterate.}
