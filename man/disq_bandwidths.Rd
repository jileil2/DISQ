\name{disq_bandwidths}
\alias{disq_bandwidths}
\title{Communication bandwidths for the manuscript design}
\usage{
disq_bandwidths(N, J, n_master, J_master, c = 0.20)
}
\arguments{
  \item{N}{Global observation count.}
  \item{J}{Global scalar constrained spline dimension.}
  \item{n_master}{Number of observations on the master worker.}
  \item{J_master}{Master worker's scalar constrained spline dimension.}
  \item{c}{Positive common multiplier, 0.20 by default.}
}
\value{A named numeric vector with the projection bandwidth \code{h}
and the master Jacobian bandwidth \code{b}.}
\description{Calculates the cube-root bandwidths used in the DISQ manuscript.
Both rules use the logarithm of two plus the global sample size.}
\details{Supply one scalar spline dimension for each mesh. The dimension is
the column count of its constrained spline basis. This calculation is
available separately from the existing \code{disq} and \code{disq_main}
workflows.}
\examples{
disq_bandwidths(N = 50000, J = 723, n_master = 1559, J_master = 31)
}
