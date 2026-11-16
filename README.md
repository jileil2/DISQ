# DISQ

DISQ fits spatial quantile regression models with constant and spatially varying coefficients over a triangulated domain.

## Installation

Install the package from source with R and a working C++ compiler.

~~~r
install.packages("remotes")
remotes::install_github("jileil2/DISQ")
~~~

DESCRIPTION lists the R dependencies and the GitHub sources for BPST, Triangulation and SQBiT. On macOS, the compiler is available through the Xcode Command Line Tools. Windows users need the Rtools version for their R release.

## Fitting a model

Use `subregion_const()` to form worker subregions from a triangulation. If the allocation of constant and varying slopes is unknown, use `forward.disq()` for forward GACV selection. Pass the chosen allocation to `disq_main()`. Its help page gives a worked example, including the spline basis and the optional evaluation grid.

~~~r
library(DISQ)
help("disq_main", package = "DISQ")
~~~

## Manuscript computations

[DISQ-analysis](https://github.com/jileil2/DISQ-analysis) is the companion repository for *Distributed Inference for Spatial Quantile Regression*. The existing package fitting functions use their original communication implementation and defaults. Use the manuscript analysis files for the paper.

`disq_bandwidths()` computes the manuscript bandwidths for a supplied design. For each mesh, give the number of columns of its constrained basis for a single coefficient function.

~~~r
disq_bandwidths(N = 50000, J = 723, n_master = 1559, J_master = 31)
~~~

### Manuscript estimator from local fits

Use `disq_from_local_fits()` to combine local Gaussian-SQBiT fits at a selected allocation. Supply the saved raw spline coefficients, constraint bases and worker rows. The function solves the local decorrelated equations, averages their estimates, applies the global smoothness correction and fits the density-weighted global projection. It holds that projection fixed during the 15 implicit master updates and evaluates the spatial covariance at the returned endpoint.

~~~r
help("disq_from_local_fits", package = "DISQ")
~~~

The function returns the actual score trace and flags incomplete updates. It takes the covariance cutoff in the supplied coordinate units. Selection and the pointwise wild bootstrap for the paper are documented with the analysis files.
