# remoteoutcome: Program Evaluation with Remotely Sensed Variables

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

R package for estimating treatment effects using remotely sensed variables (RSVs) such as satellite images or mobile phone data.

## Overview

This package implements the nonparametric methods developed in:

> Rambachan, A., Singh, R., and Viviano, D. (2025). "Program Evaluation with Remotely Sensed Outcomes." [arXiv:2411.10959](https://arxiv.org/abs/2411.10959)

**Key features:**
- Estimates treatment effects when outcomes are measured via remotely sensed variables
- Combines experimental and observational samples under a stability assumption
- Provides valid inference without rate conditions on machine learning predictions
- Supports flexible sample structures (experimental-only, observational-only, or both)
- Implements sample splitting and K-fold cross-fitting for prediction fitting

## Installation from source

```r
package_path <- "path/to/remoteoutcome"
setwd(package_path)
source("build_packages.R")
```

## Quick Start

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


## Prediction Fitting Methods

### 1. No Splitting

"Sample splitting may be eliminated under complexity restrictions that tolerate simple machine learning procedures. See e.g. Chernozhukov et al. (2020) for a recent summary."

Uses all data for both training and testing (not recommended due to overfitting).

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

### 2. Sample Splitting

Randomly splits data into training and test sets.

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

### 3. Cross-Fitting (Recommended)

K-fold cross-fitting splits data into K folds, fits predictions on K-1 folds, and predicts on the held-out fold.

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



## User-Provided Predictions

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

## Vignettes

The package includes detailed vignettes demonstrating how to use the `remoteoutcome` package to estimate treatment effects when outcomes are measured using remotely sensed variables (RSVs). Specifically, both vignettes replicate **Figure 7**  of *Rambachan, Singh, and Viviano (2025)*, with only minor differences arising from implementation updates that do **not** affect the estimator or any of its statistical properties. The updates are as follows:

1. The paper estimated \(\text{PRED}_Y(R), \text{PRED}_D(R), \text{PRED}_{S_e}(R), \text{PRED}_{S_o}(R)\) using `randomForest` package This vignette instead uses `ranger`, a modern random forest package, that improves speed and reproducibility. 

2. The paper computed the bootstrap standard errors using the `boot` package with parallelization without setting the seed. Three prediction tasks (`Ycons`, `Ylowinc`, and `Ymidinc`) were executed in parallel, each using roughly 80 cores. We now use a *custom cluster bootstrap* that supports deterministic seeding via `fix_seed = TRUE`, setting `set.seed(b)` for the *b*‑th replication. This preserves the resampling logic, ensuring reproducibility even under parallel execution.

Under these changes, both point estimates and confidence intervals remain numerically identical to those in the paper. Any remaining minute deviations stem solely from differences in random-number-generation (RNG) behavior, not from the estimator itself.

The vignettes can be accessed via:

```r
library(remoteoutcome)
vignette("treatment-effects") # A vignette that retrain the models using ranger package.
vignette("treatment-effects-original-predictions") # A vignette that uses the original randomForest-based predictions from the paper and applies the RSV estimator using those precomputed values. 
```


## Planned Features

The current version (0.1.0) implements the core RSV estimator for binary outcomes without pre-treatment covariates (Algorithm 1 from the main text). We plan to extend the package with the following features. These extensions are under active development. Contributions and feedback are welcome via [GitHub Issues](https://github.com/asheshrambachan/remoteoutcome/issues).

### 1. Discrete Outcomes

Extend the RSV estimator to handle discrete outcomes that take multiple values (not just binary 0/1). This will allow estimation of:
- Categorical treatment effects (e.g., effect on education level: none, primary, secondary, tertiary)
- Ordinal outcome variables (e.g., poverty levels: severe, moderate, mild, none)
- Count outcomes with limited support (e.g., number of assets owned: 0, 1, 2, 3+)

**Technical approach**: Generalize the moment conditions and influence functions to accommodate multi-category outcomes while maintaining the stability assumption framework.

### 2. Continuous Outcomes via Discretization

Support continuous outcome variables through automatic discretization strategies:
- Quantile-based binning (e.g., quartiles, deciles)
- Equal-width intervals
- Adaptive binning based on outcome distribution
- User-specified cutpoints

**Use cases**:
- Income or consumption expenditure (continuous)
- Test scores or health metrics (continuous)
- Any continuous outcome that can be meaningfully discretized

**Technical approach**: Pre-process continuous outcomes into discrete categories, then apply the discrete outcome RSV estimator with appropriate variance corrections for the discretization step.

### 3. Pre-treatment Covariates

Incorporate pre-treatment covariates to improve efficiency and allow for:
- Conditional average treatment effects (CATE)
- Covariate adjustment for precision gains
- Subgroup analysis
- Heterogeneous treatment effects

**Functionality**:
- `rsv_estimate(..., X = covariates)`: Add covariate matrix
- `rsv_estimate(..., X = covariates, subgroup = "gender")`: Subgroup analysis
- Doubly-robust estimation with both outcome and propensity score models

**Technical approach**: Extend Algorithm 1 to condition on pre-treatment covariates X in all conditional expectations, e.g., E[Y | R, X, S_o = 1], while maintaining the identifying stability assumption.


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
