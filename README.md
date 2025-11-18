# remoteoutcome: Program Evaluation with Remotely Sensed Variables

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

R package for estimating treatment effects using remotely sensed variables (RSVs) such as satellite images or mobile phone data.

## Overview

This package implements the nonparametric methods developed in:

> Rambachan, A., Singh, R., and Viviano, D. (2025). "Program Evaluation with Remotely Sensed Outcomes." [arXiv:2411.10959](https://arxiv.org/abs/2411.10959)

## Installation from source

The package may be installed by using the function `install_github()` from the
devtools package: 

```r
devtools::install_github(
  "https://github.com/asheshrambachan/remoteoutcome", 
  build_vignettes = TRUE
)
```

## Vignettes

The package includes two comprehensive vignettes:

### 1. Treatment Effect Estimation with Remote Sensing Variables

The package includes a vignette illustrating how the `remoteoutcome` package can be used
to estimate treatment effects when outcomes are measured using remotely sensed variables (RSVs).
In particular, it illustrates how `remoteoutcome` can be used to re-analyze the Smartcards experiment (Muralidharan et al. 2016; Muralidharan et al. 2023), following the same design reported in Section 5 of *Rambachan, Singh, and Viviano (2025)*. 

```r
vignette("treatment-effects", package = "remoteoutcome")
```

### 2. Constructing Remote Sensed Variables

This vignette provides a step-by-step guide for constructing remote sensed variables by combining multiple data sources. This vignette is essential for users who need to generate their own remote sensing features.

```r
vignette("construct-remote-vars", package = "remoteoutcome")
```

## Sample Splitting Options in `remoteoutcome`

The estimator implemented in `remoteoutcome` relies
on sample splitting. `remoteoutcome` provides the user with three possible options: 
(i) cross-fitting, (ii) sample-splitting, and (iii) no sample splitting.

### 1. Cross-Fitting (Recommended)

K-fold cross-fitting splits data into K folds, fits predictions using the 
remotely sensed variable on K-1 folds, and predicts on the held-out fold. It then
iterates across folds.

```r
library(dplyr)
library(remoteoutcome)

# Load the data
data("smartcard_data_p1", package="remoteoutcome")
data("smartcard_data_p2", package="remoteoutcome")

# Merge remote variables
smartcard_data <- inner_join(smartcard_data_p1, smartcard_data_p2, by="shrid2") %>%
rm(smartcard_data_p1, smartcard_data_p2)

data_real <- create_data_real(smartcard_data)

Y <- data_real$Ycons # binary outcome
D <- data_real$D # binary treatment
R <- data_real %>% select(starts_with("luminosity"), starts_with("satellite")) # remotely sensed variable
S_e <- !is.na(D) & (rowSums(is.na(R)) == 0) # experimental sample indicator (Observe D, R)
S_o <- !is.na(Y) & (rowSums(is.na(R)) == 0) #  observational sample indicator (Observe Y, R)
clusters <- data_real$clusters # Subdistrict-level cluster identifiers

result <- rsv_estimate(
  Y = Y, D = D, S_e = S_e, S_o = S_o, R = R,
  method = "crossfit",
  ml_params = list(nfold = 5, seed = 42),
  se_params = list(fix_seed = TRUE, clusters = clusters), 
  cores = 7
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
  ml_params = list(train_ratio = 0.5, seed = 42),
  se_params = list(fix_seed = TRUE, clusters = clusters), 
  cores = 7
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
    seed = 42             #   A random seed for each RF for reproducibility
  ),
  se = TRUE,
  se_params = list(       # Customize cluster-bootstrap standard errors:
    B = 1000,             #   Number of bootstrap replications
    clusters = clusters,  #   Cluster identifiers for clustered sampling, if not provided, use individual-level bootstrap
    fix_seed = TRUE       #   Enables deterministic seeding for reproducibility 
    ),
  cores = 7
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

# 90% confidence interval
confint(result, level = 0.90)
#>         5.0 %      95.0 %
#> D -0.03319573 0.006176239
```

## User-Provided Predictions in `remoteoutcome`

`remoteoutcome` allows the user to pass their own fitted predictions using the 
remotely sensed variable. This can be useful if the user would like to train predictions 
using the remotely sensed variable using more complex machine learning methods 
that are not directly implemented by `remoteoutcome.`

If you have your own fitted predictions, provide them directly:

```r
# Fit your own models to obtain predictions. 
# Load sample data
data("pred_real_Ycons", package = "remoteoutcome")
force(pred_real_Ycons)

result <- rsv_estimate(
  Y = pred_real_Ycons$Y,
  D = pred_real_Ycons$D,
  S_e = pred_real_Ycons$S_e,
  S_o = pred_real_Ycons$S_o,
  pred_Y = pred_real_Ycons$pred_Y,
  pred_D = pred_real_Ycons$pred_D,
  pred_S_e = pred_real_Ycons$pred_S_e,
  pred_S_o = pred_real_Ycons$pred_S_o,
  theta_init = attr(pred_real_Ycons, "theta_init"), # -0.03220447
  method = "predictions",
  se = TRUE,
  se_params = list(B = 1000, fix_seed = TRUE, clusters = pred_real_Ycons$clusters),
  cores = 7
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
#> Method: predictions
```


## Planned Features for `remoteoutcome`

The current version (0.1.0) implements the core estimator proposed in 
*Rambachan, Singh and Viviano (2025)* for binary outcomes without pre-treatment covariates 
assuming there are no direct effects. 
This corresponds to Algorithm 1 in Rambachan, Singh and Viviano (2025); see the
paper for further details. 
The following extensions are under active development: 
1. Incorporating pre-treatment covariates and discrete outcomes as discussed in Appendix E of *Rambachan, Singh and Viviano (2025)*;
2. Allowing for direct effects in the complete case (see the paper for further details); and 
3. Allowing for continuous outcomes through discretization. 

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
