# Distributed Inference for Spatial Quantile Regression (DISQ)

In this work, we develop a scalable procedure for fitting the partially linear quantile spatial model on large spatial datasets over triangulated domains.

## Quantile Spatial Model (QSM)

Spatial data often exhibit both spatial heterogeneity and tail heterogeneity. To accommodate these features, we consider the **quantile spatial model (QSM)**, in which the conditional quantile function is expressed as a combination of covariates with both constant and spatially varying effects. The varying effects are allowed to change smoothly across space.

This framework is flexible and interpretable. It allows users to study how covariate effects vary not only across spatial locations, but also across quantile levels. As a result, QSM is especially useful when the impact of predictors differs across the response distribution or when spatial patterns are not adequately captured by a purely global model.

## Distributed Inference (What it is and Why it is important?)

With the high resolution and sheer volume of modern spatial data, the computational and memory demands of statistical modeling often exceed the capacity of a single machine. This challenge arises from both the intensive linear algebra operations and the substantial storage required for statistical quantities such as the Hessian matrix. As a result, applying traditional statistical methods directly to large-scale spatial data becomes infeasible.

Distributed inference addresses this challenge by:

1. Dividing the data across multiple local machines,
2. Performing estimation locally on each machine, and
3. Aggregating the local estimates into a final result.

The primary goal of distributed inference is to achieve the same level of statistical accuracy as if the entire dataset had been processed centrally.

## Methodology (What specifically do we do?)

We propose a procedure called **Divide-and-Conquer Smoothed Quantile Bivariate Triangulation (DISQ)**, which enables the estimation of the Quantile Spatial Model (QSM) in a distributed computing environment. The method consists of three main components:

1. **Initial Estimation**  
   DISQ begins by applying the **Smoothed Quantile Bivariate Triangulation (SQBiT)** method to each local dataset. SQBiT is well suited for irregular spatial domains with moderate sample sizes. For implementation details, see the SQBiT repository.

2. **Decorrelation Step**  
   The method then orthogonalizes the constant effect estimator with respect to the varying coefficient estimator to reduce their dependency. This step ensures that the constant effect estimates are robust to the estimation error in the spatially varying components.

3. **Aggregation Step**  
   DISQ aggregates the local estimators differently for constant and varying effects.
   
   - **Constant effects**  
     DISQ performs multi-round communication between the master and local machines to refine the constant effect estimates. After a finite number of communication rounds, the aggregated estimator achieves the same convergence rate as if the entire dataset were used centrally.
   
   - **Varying effects**  
     Since varying coefficients are estimated independently in each spatial subregion, they may lack global smoothness. DISQ improves their continuity by applying a projection step across subregions.

Note that the accuracy of varying coefficient estimates cannot generally be improved via inter-machine communication because they are inherently local to each spatial partition.

## Main workflow

In practice, a typical DISQ analysis has three main steps:

1. Use `subregion_const()` to construct subregions from the global triangulation.
2. If it is unclear which covariates should have constant effects and which should have varying effects, use `forward.disq()` for model selection.
3. Fit the final distributed model using `disq_main()`.

## Inference

DISQ supports inference for both constant and varying coefficients.

1. **Constant coefficients**  
   DISQ provides confidence intervals through the communication-based aggregation step.

2. **Varying coefficients**  
   DISQ optionally supports distributed wild bootstrap inference for varying coefficients through the `vc.wb = TRUE` option in `disq_main()`.

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
   Fits the final DISQ model by combining local SQBiT estimation, decorrelation-based aggregation for constant coefficients, and smoothness correction for varying coefficients.

### How to install the packages in R?

```r
install.packages(c("Rcpp", "RcppArmadillo", "devtools"))

devtools::install_github("FIRST-Data-Lab/BPST", force = TRUE)
devtools::install_github("funstatpackages/Triangulation", force = TRUE)
devtools::install_github("jileil2/SQBiT", force = TRUE)
devtools::install_github("jileil2/DISQ", force = TRUE)