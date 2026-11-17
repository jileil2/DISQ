# DISQ: Distributed Inference for Spatial Quantile Regression

`DISQ` fits spatial quantile regression models with constant and spatially
varying coefficients. Local bivariate spline fits describe variation over a
triangulated domain. Decorrelation and iterative communication combine
information across subregions for inference on the constant coefficients.

The package contains domain decomposition, local estimation, forward GACV
selection, communication, and spatial covariance routines. The companion
[DISQ-analysis](https://github.com/jileil2/DISQ-analysis) repository contains
the manuscript simulations, CHD analysis, and figure and table scripts.

## Installation

Use R 3.6 or newer and a C++ compiler. On macOS, install the Xcode Command Line
Tools; on Windows, use the Rtools version for your R release.

```r
install.packages("remotes")
remotes::install_github("jileil2/DISQ")
library(DISQ)
```

The installer obtains SQBiT, BPST, and Triangulation from the GitHub sources
listed in [DESCRIPTION](DESCRIPTION). To install a local checkout, run
`R CMD INSTALL .` from its root after installing the dependencies.

## Getting started

The [worked example](inst/examples/getting-started.R) simulates a spatial
response, fits local SQBiT models on two subregions, and combines them using
a fixed global decorrelating projection. It runs 15 communication rounds and
returns constant-coefficient estimates, spatial intervals, and fitted quantiles.
Run it in R after installation:

```r
source(system.file("examples", "getting-started.R", package = "DISQ"))
```

For an analysis with saved local fits, supply one entry per worker in `local`:

```r
fit <- disq_from_local_fits(
  Y = Y, C = C, X = X, S = S,
  B = B, Q2 = Q2, index = index, local = local,
  tau = 0.5, radius = 0.035, c = 0.20, rounds = 15, level = 0.01
)

fit$eta                       # Constant coefficients.
fit$covariance$ci             # 99% spatial covariance intervals.
fit$beta                      # Smoothness-corrected varying coefficients.
fit$fitted                    # Fitted conditional quantiles.
fit$communication$trace       # Initialization and each communication update.
```

`C` contains the constant terms; `X` contains the varying terms, including a
column of ones for a varying intercept. `B` is the raw global spline basis,
`Q2` imposes its smoothness constraints, and `index` gives the observation rows
for each worker. Each local entry contains its raw spline column indices,
constraint basis, spline coefficients, constant estimates, and local bandwidth.
See `?disq_from_local_fits` for the input format and returned diagnostics.
The covariance radius uses the units of the supplied spatial coordinates.

## Analysis workflow

1. Define the triangulation and build its spline basis with BPST.
2. Use `subregion_const()` to group triangles into worker subregions.
3. If the coefficient allocation is unknown, use `forward.disq()` to compare
   constant and varying slopes by distributed GACV.
4. Fit local models at the chosen allocation and combine the results with
   `disq_from_local_fits()`.

`disq_bandwidths()` gives the global projection and master bandwidths. Its
default common multiplier is `c = 0.20`; `J` and `J_master` are the constrained
spline dimensions for a single coefficient function.

```r
disq_bandwidths(N = 50000, J = 723, n_master = 1559, J_master = 31)
```

The earlier `disq_main()` interface fits local models, calls `disq()`, and
optionally runs a varying-coefficient wild bootstrap. It retains its original
communication routine and defaults. Use `disq_from_local_fits()` for the
fixed-projection communication estimator in the current manuscript, and the
companion analysis repository for the manuscript's selection and bootstrap
calculations.

## Package files

| Directory | Contents |
|---|---|
| `R/` | Fitting, selection, partitioning, projection, communication, and inference |
| `src/` | C++ matrix routines and Rcpp bindings |
| `man/` | R help pages |
| `inst/examples/` | A self-contained analysis example |
| `tests/` | Numerical checks for bandwidths, projection, communication, and covariance |

Run `R CMD build .` and `R CMD check --no-manual DISQ_0.2.1.tar.gz` to check a
source release. [NEWS.md](NEWS.md) records package changes.

## License

MIT. See [LICENSE](LICENSE).
