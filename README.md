# Distributed Inference for Spatial Quantile Regression (DISQ)

In this work, we develop a scalable procedure for fitting the partially linear quantile spatial model on large spatial datasets over triangulated domains.

## Quantile Spatial Model (QSM)

Spatial data often exhibit both spatial heterogeneity and tail heterogeneity. To accommodate these features, we consider the **quantile spatial model (QSM)**, in which the conditional quantile function is expressed as a combination of covariates with both constant and spatially varying effects. The varying effects are allowed to change smoothly across space.

This framework is flexible and interpretable. It allows users to study how covariate effects vary not only across spatial locations, but also across quantile levels. As a result, QSM is especially useful when the impact of predictors differs across the response distribution or when spatial patterns are not adequately captured by a purely global model.

## Methodology

We propose **Distributed Inference for Spatial Quantile Regression (DISQ)**, a divide-and-combine procedure for fitting QSM on large spatial datasets.

DISQ combines three main ideas:

1. **Domain decomposition.** The global triangulated domain is partitioned into smaller subregions.
2. **Local estimation.** On each subregion, the model is fitted locally using the SQBiT method.
3. **Aggregation.** The local estimators are combined to recover global estimators of both constant and varying coefficients.

This strategy offers several advantages.

First, DISQ is computationally scalable for large datasets, since the full problem is broken into smaller subproblems that can be solved in parallel.

Second, the method remains well suited for irregular spatial domains because it is built on bivariate triangulation splines.

Third, DISQ treats constant and varying coefficients differently in the aggregation step. The varying coefficients are combined through a smoothness correction so that the final estimator is globally smooth across subregion boundaries. The constant coefficients are further refined through a communication-based aggregation algorithm to improve statistical efficiency.

## Main workflow

In practice, a typical DISQ analysis has three main steps:

1. Use `subregion_const()` to construct subregions from the global triangulation.
2. If it is unclear which covariates should have constant effects and which should have varying effects, use `forward.disq()` for model selection.
3. Fit the final distributed model using `disq_main()`.

## Inference

DISQ supports inference for both constant and varying coefficients.

1. **Constant coefficients.** DISQ provides confidence intervals through the communication-based aggregation step.
2. **Varying coefficients.** DISQ optionally supports distributed wild bootstrap inference for varying coefficients through the `vc.wb = TRUE` option in `disq_main()`.

## Instructions for use

### Dependencies

Our codes depend heavily on the packages **Rcpp**, **RcppArmadillo**, **BPST**, **Triangulation**, and **SQBiT**.

In particular:

- **Rcpp** and **RcppArmadillo** are used for efficient linear algebra and computational routines.
- **BPST** is used for bivariate penalized splines over triangulations.
- **Triangulation** is used to construct triangulation meshes over spatial domains.
- **SQBiT** is used for local estimation within each subregion.

The packages **Rcpp** and **RcppArmadillo** can be installed from CRAN using `install.packages()`. The packages **BPST**, **Triangulation**, and **SQBiT** can be installed from GitHub.

### Main functions

The main functions in this package are:

1. `subregion_const`  
   Constructs subregions from the global triangulation for distributed estimation.

2. `forward.disq`  
   Performs forward model selection within the DISQ framework. This function is recommended when users do not know which covariates should receive constant coefficients and which should receive varying coefficients.

3. `disq_main`  
   Fits the final DISQ model by combining local SQBiT estimation, smoothness correction for varying coefficients, and communication-based aggregation for constant coefficients.

### How to install the packages in R?

```r
install.packages(c("Rcpp", "RcppArmadillo", "devtools"))

devtools::install_github("FIRST-Data-Lab/BPST", force = TRUE)
devtools::install_github("funstatpackages/Triangulation", force = TRUE)
devtools::install_github("jileil2/SQBiT", force = TRUE)