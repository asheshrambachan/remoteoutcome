# remoteoutcome

**Program Evaluation with Remotely Sensed Outcomes**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`remoteoutcome` implements the RSV (Remotely Sensed Variable) estimator of
Rambachan, Singh, and Viviano (2025) for estimating average treatment effects
when outcomes are imperfectly measured by a remotely sensed variable — such as
satellite imagery, nighttime luminosity, or mobile phone activity.

The estimator combines an experimental sample (treatment assigned, outcome
possibly unobserved) with an observational sample (outcome observed, treatment
not assigned) by using the remote sensing variable as a bridge between the two
sources of variation. It supports binary and multi-category outcomes,
cross-fitting, and both influence-function and score-bootstrap standard errors.

## Installation

Install from CRAN:

```r
install.packages("remoteoutcome")
```

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("asheshrambachan/remoteoutcome", build_vignettes = TRUE)
```

## Quick start

```r
library(remoteoutcome)

# Simulate data: 300 experimental, 700 observational, 100 overlap units
dat <- sim_rsv_data(n_e = 300, n_o = 700, n_v = 100, tau = 0.10, seed = 42)
R   <- as.matrix(dat[, paste0("R", 1:5)])

# Fit RSV estimator with 2-fold cross-fitting (logistic nuisance models)
fit <- cv.rsv(
  Y        = dat$Y,
  D        = dat$D,
  S_e      = dat$S_e,
  S_o      = dat$S_o,
  R        = R,
  y_levels = c(0L, 1L),
  models   = list(
    Y   = list(model = "logit"),
    D   = list(model = "logit"),
    S_e = list(model = "logit"),
    S_o = list(model = "logit")
  ),
  nfolds = 2,
  seed   = 42L
)

# Add standard errors (analytical influence function)
fit <- add_se(fit, method = "influence")
summary(fit)
```

## Vignettes

- **Getting started** — introduces the estimator, data structure, and basic workflow
- **Replication: Smartcards in India** — replicates the empirical analysis from Rambachan, Singh, and Viviano (2025) using the Muralidharan et al. (2016) smartcards data
- **Replication: Crop Burning in India** — replicates the empirical analysis from Rambachan, Singh, and Viviano (2025) using the Jack et al. (2025) crop burning data
- **Replication: Uganda Forest Cover** — replicates the simulation study from Rambachan, Singh, and Viviano (2025) calibrated to Jayachandran et al. (2017)

```r
browseVignettes("remoteoutcome")
```

## Citation

If you use this package, please cite:

```bibtex
@article{rambachan2025program,
  title   = {Program Evaluation with Remotely Sensed Outcomes},
  author  = {Rambachan, Ashesh and Singh, Rahul and Viviano, Davide},
  journal = {arXiv preprint arXiv:2411.10959},
  year    = {2025}
}
```

## Issues and contributions

Please report issues or suggest improvements at the
[GitHub repository](https://github.com/asheshrambachan/remoteoutcome/issues).

## License

MIT — see [LICENSE](LICENSE) for details.
