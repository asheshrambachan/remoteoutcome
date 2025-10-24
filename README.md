# remoteoutcome: Program Evaluation with Remotely Sensed Variables

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

R package for estimating treatment effects using remotely sensed variables (RSVs) such as satellite images or mobile phone data.

## Overview

This package implements the nonparametric methods developed in:

> Rambachan, A., Singh, R., and Viviano, D. (2025). "Program Evaluation with Remotely Sensed Outcomes." [arXiv:2411.10959](https://arxiv.org/abs/2411.10959)

## Installation from source

The package may be installed by using the function `install_github()` from the
remotes package: 

```r
# Install remotes package if not installed
install.packages("remotes")

# Turn off warning-error-conversion, because the tiniest warning stops installation
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

# install from github
remotes::install_github("https://github.com/asheshrambachan/remoteoutcome")
```

## Vignettes

The package includes vignettes demonstrating how to use the `remoteoutcome` package 
to estimate treatment effects when outcomes are measured using remotely sensed variables (RSVs). 
The vignettes replicate the re-analysis of the Smartcards experiment (Muralidharan et al. 2016; Muralidharan et al. 2023) reported in Section 5 of *Rambachan, Singh, and Viviano (2025)* (specifically, Figure 7). 
The vignettes replicate the analysis in two ways:

1. `treatment-effects-original-predictions.Rmd:` uses the precomputed random forest
predictions based on the remotely sensed variables constructed in *Rambachan, 
Singh, and Viviano (2025).*
2. `treatment-effects.Rmd:` retrains the random forest predictions based on the 
remotely sensed variables. 

The vignettes can be accessed via:

```r
library(remoteoutcome)
vignette("treatment-effects") # A vignette that retrain the models using ranger package.
vignette("treatment-effects-original-predictions") # A vignette that uses the original randomForest-based predictions from the paper and applies the RSV estimator using those precomputed values. 
```

## Sample Splitting Options in `remoteoutcome`

The estimator for the average treatment effect in the experimental sample relies
on sample splitting. `remoteoutcome` provides the user with three possible options: 
(i) cross-fitting, (ii) sample-splitting, and (iii) no sample splitting.

### 1. Cross-Fitting (Recommended)

K-fold cross-fitting splits data into K folds, fits predictions using the 
remotely sensed variable on K-1 folds, and predicts on the held-out fold. It then
iterates across folds.

```r
result <- rsv_estimate(
  Y = Y, D = D, S_e = S_e, S_o = S_o, R = R,
  method = "crossfit",
  ml_params = list(nfold = 5, seed = 42, cores = 7),
  se_params = list(fix_seed = TRUE, clusters = clusters)
)

print(result)
#> RSV Treatment Effect Estimate
#> ==============================
#> 
#> Coefficient: -0.1082
#> 
#> Method: crossfit
```

### 2. Sample Splitting

Sample splitting splits the data into a train and test set. The train set is used
to fit predictions using the remotely sensed variable, and the test set is used
to construct the estimator. 

```r
result <- rsv_estimate(
  Y = Y, D = D, S_e = S_e, S_o = S_o, R = R,
  method = "split",
  ml_params = list(train_ratio = 0.5, seed = 42, cores = 7),
  se_params = list(fix_seed = TRUE, clusters = clusters)
)

print(result)
#> RSV Treatment Effect Estimate
#> ==============================
#> 
#> Coefficient: -0.1086 (SE: 0.0959)
#> 
#> Sample sizes:
#>   Experimental: 3032
#>   Observational: 2575
#>   Both: 1451
#> 
#> Method: split
```

### 3. No Splitting

No sample splitting uses all of the data for both training predictions based 
on the remotely sensed variable and constructing the estimator. This is 
generally not recommended due to possible overfitting, but valid estimation/inference
can still be conducted without sample splitting for particularly simple 
machine learning procedures. See Rambachan, Singh and Viviano (2025) for more discussion.

```r
result <- rsv_estimate(
  Y = Y, D = D, S_e = S_e, S_o = S_o, R = R,
  method = "none",
  ml_params = list(seed = 42, cores = 7),
  se_params = list(fix_seed = TRUE, clusters = clusters, cores = 7)
)

print(result)
#> RSV Treatment Effect Estimate
#> ==============================
#> 
#> Coefficient: -0.0135 (SE: 0.0120)
#> 
#> Sample sizes:
#>   Experimental: 6055
#>   Observational: 5186
#>   Both: 2929
#> 
#> Method: none
```

## User-Provided Predictions in `remoteoutcome`

`remoteoutcome` allows the user to pass their own fitted predictions using the 
remotely sensed variable. This can be useful if the user would like to train predictions 
using the remotely sensed variable using more complex machine learning methods 
that are not directly implemented by `remoteoutcome.`

If you have your own fitted predictions, provide them directly:

```r
# Fit your own models to obtain predictions. 
# Example data: pred_real_Ycons (included in package)
result <- rsv_estimate(
  Y = pred_real_Ycons$Y,
  D = pred_real_Ycons$D,
  S_e = pred_real_Ycons$S_e,
  S_o = pred_real_Ycons$S_o,
  pred_Y = pred_real_Ycons$pred_Y,
  pred_D = pred_real_Ycons$pred_D,
  pred_S_e = pred_real_Ycons$pred_S_e,
  pred_S_o = pred_real_Ycons$pred_S_o,
  se = TRUE,
  se_params = list(B = 1000, fix_seed = TRUE, cores = 7, clusters = pred_real_Ycons$clusters),
)

print(result)
#> RSV Treatment Effect Estimate
#> ==============================
#> 
#> Coefficient: -0.0127 (SE: 0.0095)
#> 
#> Sample sizes:
#>   Experimental: 6055
#>   Observational: 5186
#>   Both: 2929
#> 
#> Method: none
```

## Quick Start Example of `remoteoutcome`

We provide a brief illustration of how ``remoteoutcome'' can be used to estimate
treatment effects. 

#### Load sample data

```r
library(dplyr)
library(remoteoutcome)

# Example data: data_real (included in package)
Y <- data_real$Ycons # binary outcome
D <- data_real$D # binary treatment
R <- data_real %>% select(starts_with("luminosity"), starts_with("satellite")) # remotely sensed variable
S_e <- !is.na(D) & (rowSums(is.na(R)) == 0) # experimental sample indicator (Observe D, R)
S_o <- !is.na(Y) & (rowSums(is.na(R)) == 0) #  observational sample indicator (Observe Y, R)
clusters <- data_real$clusters # Subdistrict-level cluster identifiers
```

### Basic Example

```r
# Estimate treatment effect
result <- rsv_estimate(
  Y = Y,
  D = D,
  S_e = S_e,
  S_o = S_o,
  R = R,
  eps = 1e-2,
  method = "none",
  ml_params = list(       # Customize random forest parameters:
    ntree = 100,          #   Number of trees
    classwt_Y = c(10, 1), #   Class weights for PRED_Y model
    seed = 42,            #   A random seed for each RF for reproducibility
    cores = 7             #   Number of cores used for training
  ),
  se = TRUE,
  se_params = list(       # Customize cluster-bootstrap standard errors:
    B = 1000,             #   Number of bootstrap replications
    clusters = clusters,  #   Cluster identifiers for clustered sampling, if not provided, use individual-level bootstrap
    fix_seed = TRUE,      #   Enables deterministic seeding for reproducibility 
    cores = 7             #   Number of cores for bootstrap replications
    )
)

summary(result)
#> RSV Treatment Effect Estimate
#> ==============================
#> 
#> Coefficient:
#>   Estimate Std.Error t.value Pr..t.
#> D -0.01351   0.01197  -1.129  0.259
#> 
#> Sample sizes:
#>   Experimental: 6055
#>   Observational: 5186
#>   Both: 2929
#> 
#> Prediction fitting method: none
#> 
#> Call:
#> rsv_estimate(Y = Y, D = D, S_e = S_e, S_o = S_o, R = R, eps = 0.01, 
#>     method = "none", ml_params = list(ntree = 100, classwt_Y = c(10, 
#>         1), seed = 42, cores = 7), se = TRUE, se_params = list(B = 1000, 
#>         clusters = clusters, fix_seed = TRUE, cores = 7))

# 90% confidence interval
confint(result, level = 0.90)
#>         5.0 %      95.0 %
#> D -0.03319573 0.006176239
```

## Planned Features for `remoteoutcome`

The current version (0.1.0) implements the core estimator proposed in 
Rambachan, Singh and Viviano (2025) for binary outcomes without pre-treatment covariates 
assuming there are no direct effects. 
This corresponds to Algorithm 1 in Rambachan, Singh and Viviano (2025); see the
paper for further details. 
The following extensions are under active development: (1) incorporating
pre-treatment covariates and discrete outcomes as discussed in Appendix E of Rambachan, Singh and Viviano (2025);
(2) allowing for direct effects in the complete case (see the paper for further details); and 
(3) allowing for continuous outcomes through discretization. 

Contributions and feedback are welcome via [GitHub Issues](https://github.com/asheshrambachan/remoteoutcome/issues).

## Citation

If you use this package, please cite:

```bibtex
@article{rambachan2025program,
  title={Program Evaluation with Remotely Sensed Outcomes},
  author={Rambachan, Ashesh and Singh, Rahul and Viviano, Davide},
  journal={arXiv preprint arXiv:2411.10959},
  year={2025}
}
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Issues and Contributions

Please report issues or suggest improvements at the [GitHub repository](https://github.com/asheshrambachan/remoteoutcome/issues).
