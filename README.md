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

## Installation

```r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install remoteoutcome package from source
devtools::install_local("path/to/remoteoutcome")
```

### Dependencies

The package requires:
- R >= 4.0.0
- dplyr >= 1.0.0
- randomForest >= 4.6-14
- boot >= 1.3-28
- fixest >= 0.10.0

## Quick Start

### Basic Example

```r
library(remoteoutcome)

# Example data structure:
# Y: binary outcome (NA for experimental-only observations)
# D: binary treatment (NA for observational-only observations)
# S: sample indicator ("e" = experimental, "o" = observational, "both" = both)
# R: remotely sensed variable (e.g., satellite image features)

# Estimate treatment effect with 5-fold cross-fitting
result <- rsv_estimate(
  Y = Y,
  D = D,
  S = S,
  R = R,
  method = "crossfit",
  nfolds = 5,
  se = TRUE,
  alpha = 0.05  # 95% confidence interval
)

# View results
print(result)
#> RSV Treatment Effect Estimate
#> ==============================
#>
#> Coefficient: 0.1234 (SE: 0.0456)
#> 95% CI: [0.0340, 0.2128]
#>
#> Sample sizes:
#>   Experimental: 1000
#>   Observational: 800

summary(result)

# Extract coefficient
coef(result)

# Extract confidence interval
confint(result)
```

## Sample Indicator Options

The package handles three types of sample membership (Remark 3 in the paper):

| Sample Type | Code | Observations |
|------------|------|--------------|
| Experimental only | `S = "e"` | Observe D, R (not Y) |
| Observational only | `S = "o"` | Observe Y, R (not D) |
| Both samples | `S = "both"` | Observe Y, D, R |

**Example:**
```r
# Create sample indicator
S <- rep("e", 1000)  # Experimental sample
S[501:800] <- "o"    # Observational sample
S[801:1000] <- "both"  # In both samples

# Y is NA for experimental-only observations
Y[S == "e"] <- NA

# D is NA for observational-only observations
D[S == "o"] <- NA
```

## Prediction Fitting Methods

### 1. Cross-Fitting (Recommended)

K-fold cross-fitting splits data into K folds, fits predictions on K-1 folds, and predicts on the held-out fold.

```r
result <- rsv_estimate(
  Y = Y, D = D, S = S, R = R,
  method = "crossfit",
  nfolds = 5,
  se = TRUE
)
```

### 2. Sample Splitting

Randomly splits data into training and test sets.

```r
result <- rsv_estimate(
  Y = Y, D = D, S = S, R = R,
  method = "split",
  split_ratio = 0.5,  # 50% train, 50% test
  se = TRUE
)
```

### 3. No Splitting

Uses all data for both training and testing (not recommended due to overfitting).

```r
result <- rsv_estimate(
  Y = Y, D = D, S = S, R = R,
  method = "none",
  se = TRUE
)
```

## Customization

### Custom Machine Learning Parameters

Customize random forest parameters via `ml_params`:

```r
result <- rsv_estimate(
  Y = Y, D = D, S = S, R = R,
  method = "crossfit",
  ml_params = list(
    ntree = 1000,        # Number of trees (default: 500)
    classwt_Y = c(5, 1)  # Class weights for outcome model (default: c(10, 1))
  ),
  se = TRUE
)
```

### Custom Confidence Levels

Specify confidence level via `alpha` parameter:

```r
# 90% confidence interval
result <- rsv_estimate(Y = Y, D = D, S = S, R = R, alpha = 0.10)

# 95% confidence interval (default)
result <- rsv_estimate(Y = Y, D = D, S = S, R = R, alpha = 0.05)

# 99% confidence interval
result <- rsv_estimate(Y = Y, D = D, S = S, R = R, alpha = 0.01)
```

### Clustered Standard Errors

For clustered data, provide cluster identifiers:

```r
result <- rsv_estimate(
  Y = Y, D = D, S = S, R = R,
  method = "crossfit",
  clusters = cluster_id,  # Vector of cluster identifiers
  B = 1000,              # Number of bootstrap replications
  se = TRUE
)
```

## User-Provided Predictions

If you have your own fitted predictions, provide them directly:

```r
# Fit your own models to obtain:
# pred_Y: E[Y | R, S_o = 1]
# pred_D: E[D | R, S_e = 1]
# pred_S_e: P(S_e = 1 | R)
# pred_S_o: P(S_o = 1 | R)

result <- rsv_estimate(
  Y = Y,
  D = D,
  S = S,
  pred_Y = pred_Y,
  pred_D = pred_D,
  pred_S_e = pred_S_e,
  pred_S_o = pred_S_o,
  se = TRUE
)
```

## Reproducibility

Set a random seed for reproducibility in sample splitting/cross-fitting:

```r
result <- rsv_estimate(
  Y = Y, D = D, S = S, R = R,
  method = "crossfit",
  seed = 12345,
  se = TRUE
)
```

## Complete Example

```r
library(remoteoutcome)

# Generate synthetic data
set.seed(123)
n <- 1000

# Remotely sensed variable
R <- rnorm(n)

# True outcome (latent)
Y_latent <- rbinom(n, 1, plogis(0.5 * R))

# Treatment
D <- rbinom(n, 1, plogis(-0.3 * R))

# Sample indicator
S <- sample(c("e", "o", "both"), n, replace = TRUE, prob = c(0.4, 0.4, 0.2))

# Observed outcome (NA for experimental-only)
Y <- Y_latent
Y[S == "e"] <- NA

# Observed treatment (NA for observational-only)
D_obs <- D
D_obs[S == "o"] <- NA

# Estimate treatment effect
result <- rsv_estimate(
  Y = Y,
  D = D_obs,
  S = S,
  R = R,
  method = "crossfit",
  nfolds = 5,
  ml_params = list(ntree = 500, classwt_Y = c(10, 1)),
  se = TRUE,
  clusters = NULL,
  B = 1000,
  alpha = 0.05,
  seed = 456
)

# View results
print(result)
summary(result)
coef(result)
confint(result)
```

## Functions

### Main Functions

- `rsv_estimate()`: Main estimation function
- `cv_rsv()`: Cross-validation wrapper
- `generate_rsv_data()`: Generate synthetic data for testing

### S3 Methods

- `print.rsv()`: Print results
- `summary.rsv()`: Detailed summary with inference
- `coef.rsv()`: Extract coefficient
- `vcov.rsv()`: Extract variance-covariance matrix
- `confint.rsv()`: Extract confidence interval

## Vignettes

The package includes detailed vignettes demonstrating practical applications:

- **Poverty Application** (`vignette("poverty-application", package = "remoteoutcome")`): Replicates the poverty analysis from Rambachan, Singh, and Viviano (2025) using data from an anti-poverty program in India. Shows how to:
  - Estimate treatment effects with real and synthetic sample definitions
  - Use satellite imagery features as remotely sensed variables
  - Compare RSV estimates with traditional benchmark approaches
  - Create publication-quality visualizations

## Planned Features

The current version (0.1.0) implements the core RSV estimator for binary outcomes without pre-treatment covariates (Algorithm 1 from the main text). We plan to extend the package with the following features:

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

### Implementation Timeline

These extensions are under active development. Contributions and feedback are welcome via [GitHub Issues](https://github.com/asheshrambachan/remoteoutcome/issues).

**Priority order**:
1. Pre-treatment covariates (highest priority - most requested feature)
2. Continuous outcomes via discretization
3. General discrete outcomes

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

## Authors

- **Ashesh Rambachan** - Stanford University
- **Rahul Singh**
- **Davide Viviano**

## Issues and Contributions

Please report issues or suggest improvements at the [GitHub repository](https://github.com/asheshrambachan/remoteoutcome/issues).
