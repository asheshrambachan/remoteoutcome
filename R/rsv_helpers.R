# =============================================================================
# RSV Helper Functions (Binary Outcome, No Covariates)
# Following main text Algorithm 1 with Remark 3 generalization
# =============================================================================

#' Convert sample indicator to binary indicators
#'
#' Converts sample indicator S to binary indicators for experimental and
#' observational samples. Handles the general case from Remark 3 where
#' observations can be in experimental only, observational only, or both.
#'
#' @param S Sample indicator. Can be:
#'   \describe{
#'     \item{"e"}{Experimental sample only (observe D, R)}
#'     \item{"o"}{Observational sample only (observe Y, R)}
#'     \item{"both"}{Both samples (observe Y, D, R)}
#'   }
#'
#' @return A list with:
#'   \describe{
#'     \item{S_e}{Binary indicator for experimental sample (D available)}
#'     \item{S_o}{Binary indicator for observational sample (Y available)}
#'   }
#'
#' @keywords internal
parse_sample_indicator <- function(S) {
  S_e <- (S == "e") | (S == "both")
  S_o <- (S == "o") | (S == "both")

  list(S_e = as.numeric(S_e), S_o = as.numeric(S_o))
}


#' Fit predictions using random forest
#'
#' Fits random forest models to predict Y(R), D(R), S_e(R), and S_o(R) using
#' the remotely sensed variable R as the only predictor.
#'
#' @param R Remotely sensed variable (continuous or discrete).
#' @param Y Outcome variable (binary, NA for experimental-only sample).
#' @param D Treatment indicator (binary, NA for observational-only sample).
#' @param S Sample indicator ("e", "o", or "both").
#' @param R_pred Vector of R values on which to make predictions (defaults to R).
#' @param ml_params List of parameters for random forest. Should include:
#'   \describe{
#'     \item{ntree}{Number of trees}
#'     \item{classwt_Y}{Class weights for outcome model}
#'   }
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{pred_Y}{E[Y | R, S_o = 1]}
#'     \item{pred_D}{E[D | R, S_e = 1]}
#'     \item{pred_S_e}{P(S_e = 1 | R)}
#'     \item{pred_S_o}{P(S_o = 1 | R)}
#'   }
#'
#' @keywords internal
fit_predictions_rf <- function(R, Y, D, S, R_pred = NULL,
                               ml_params = list(ntree = 500, classwt_Y = c(10, 1))) {
  if (is.null(R_pred)) R_pred <- R

  # Extract ML parameters
  ntree <- ml_params$ntree
  classwt_Y <- ml_params$classwt_Y

  # Convert R to matrix if it's a vector
  R_mat <- as.matrix(R)
  R_pred_mat <- as.matrix(R_pred)

  # Parse sample indicators
  S_parsed <- parse_sample_indicator(S)
  S_e <- S_parsed$S_e
  S_o <- S_parsed$S_o

  # Fit E[Y | R, S_o = 1] on observations with Y available
  obs_idx <- (S_o == 1)
  if (sum(obs_idx) == 0) stop("No observations in observational sample")

  model_Y <- randomForest::randomForest(
    x = R_mat[obs_idx, , drop = FALSE],
    y = factor(Y[obs_idx], levels = c(0, 1)),
    classwt = classwt_Y,
    ntree = ntree
  )

  # Fit E[D | R, S_e = 1] on observations with D available
  exp_idx <- (S_e == 1)
  if (sum(exp_idx) == 0) stop("No observations in experimental sample")

  model_D <- randomForest::randomForest(
    x = R_mat[exp_idx, , drop = FALSE],
    y = factor(D[exp_idx], levels = c(0, 1)),
    ntree = ntree
  )

  # Fit P(S_e = 1 | R) on full sample
  model_S_e <- randomForest::randomForest(
    x = R_mat,
    y = factor(S_e, levels = c(0, 1)),
    ntree = ntree
  )

  # Fit P(S_o = 1 | R) on full sample
  model_S_o <- randomForest::randomForest(
    x = R_mat,
    y = factor(S_o, levels = c(0, 1)),
    ntree = ntree
  )

  # Generate predictions
  data.frame(
    pred_Y = predict(model_Y, R_pred_mat, type = "prob")[, "1"],
    pred_D = predict(model_D, R_pred_mat, type = "prob")[, "1"],
    pred_S_e = predict(model_S_e, R_pred_mat, type = "prob")[, "1"],
    pred_S_o = predict(model_S_o, R_pred_mat, type = "prob")[, "1"]
  )
}


#' Compute RSV estimator from predictions
#'
#' Computes the RSV treatment effect estimator given predictions and sample
#' indicators. Implements Algorithm 1 from main text with Remark 3 generalization
#' allowing observations to be in experimental only, observational only, or both.
#'
#' @param Y Outcome variable (binary, NA where not observed).
#' @param D Treatment indicator (binary, NA where not observed).
#' @param S Sample indicator ("e", "o", or "both").
#' @param pred_Y Predicted E[Y | R, S_o = 1].
#' @param pred_D Predicted E[D | R, S_e = 1].
#' @param pred_S_e Predicted P(S_e = 1 | R).
#' @param pred_S_o Predicted P(S_o = 1 | R).
#' @param eps Small constant for numerical stability. Default is 1e-6.
#'
#' @return A list containing:
#'   \describe{
#'     \item{coef}{Treatment effect estimate}
#'     \item{weights}{Efficient weights used in estimation}
#'     \item{n_obs}{Sample size in observational sample}
#'     \item{n_exp}{Sample size in experimental sample}
#'     \item{n_both}{Sample size in both samples}
#'   }
#'
#' @keywords internal
rsv_compute <- function(Y, D, S, pred_Y, pred_D, pred_S_e, pred_S_o, eps = 1e-6) {

  # Parse sample indicators
  S_parsed <- parse_sample_indicator(S)
  S_e <- S_parsed$S_e  # 1 if D observed
  S_o <- S_parsed$S_o  # 1 if Y observed

  # Sample sizes
  n_exp <- sum(S_e == 1)
  n_obs <- sum(S_o == 1)
  n_both <- sum(S_e == 1 & S_o == 1)
  n <- length(S)

  # Propensity scores (bound away from 0)
  pi_e <- pmax(pred_S_e, eps)
  pi_o <- pmax(pred_S_o, eps)

  # Bound predictions away from 0 and 1
  pred_Y <- pmax(pmin(pred_Y, 1 - eps), eps)
  pred_D <- pmax(pmin(pred_D, 1 - eps), eps)

  # Compute moment functions
  # For experimental sample: mu_D(R) = E[D | R, S_e = 1]
  mu_D <- pred_D

  # For observational sample: mu_Y(R) = E[Y | R, S_o = 1]
  mu_Y <- pred_Y

  # Compute influence functions (centered moment functions)
  # psi_D = S_e * (D - mu_D) / pi_e
  # psi_Y = S_o * (Y - mu_Y) / pi_o
  psi_D <- S_e * (D - mu_D) / pi_e
  psi_Y <- S_o * (Y - mu_Y) / pi_o

  # Replace NA values with 0 (for observations where Y or D not observed)
  psi_D[is.na(psi_D)] <- 0
  psi_Y[is.na(psi_Y)] <- 0

  # Compute variance terms for efficient weighting
  # Var(D | R, S_e = 1) and Var(Y | R, S_o = 1)
  var_D <- mu_D * (1 - mu_D)
  var_Y <- mu_Y * (1 - mu_Y)

  # Initial estimate for theta (needed for weight construction)
  # Use simple difference in means as initial value
  theta_init <- mean(D[S_e == 1], na.rm = TRUE) -
                mean(D[S_e == 1 & D == 0], na.rm = TRUE)
  if (is.na(theta_init)) theta_init <- 0

  # Compute sigma^2 for optimal weighting (Equation in main text)
  # sigma^2(R) = Var(D | R) / pi_e^2 + theta^2 * Var(Y | R) / pi_o^2
  sigma2 <- var_D / pi_e^2 + theta_init^2 * var_Y / pi_o^2
  sigma2 <- pmax(sigma2, eps)

  # Efficient weight function
  # h(R) = Var(Y | R) / (pi_o^2 * sigma^2(R))
  h <- (var_Y / pi_o^2) / sigma2

  # Compute weighted estimator
  # theta = E[psi_D * h] / E[psi_Y * h]
  numerator <- mean(psi_D * h)
  denominator <- mean(psi_Y * h)

  coef <- numerator / denominator

  list(
    coef = coef,
    weights = h,
    n_obs = n_obs,
    n_exp = n_exp,
    n_both = n_both,
    numerator = numerator,
    denominator = denominator
  )
}


#' Bootstrap standard errors for RSV estimator
#'
#' Computes cluster-bootstrap standard errors for the RSV estimator.
#'
#' @param Y Outcome variable.
#' @param D Treatment indicator.
#' @param S Sample indicator.
#' @param pred_Y Predicted E[Y | R, S_o = 1].
#' @param pred_D Predicted E[D | R, S_e = 1].
#' @param pred_S_e Predicted P(S_e = 1 | R).
#' @param pred_S_o Predicted P(S_o = 1 | R).
#' @param clusters Cluster identifiers for bootstrap.
#' @param B Number of bootstrap replications. Default is 1000.
#' @param alpha Significance level for confidence interval. Default is 0.05 (95% CI).
#' @param eps Small constant for numerical stability.
#'
#' @return A list containing:
#'   \describe{
#'     \item{se}{Standard error of treatment effect}
#'     \item{ci_lower}{Lower bound of confidence interval}
#'     \item{ci_upper}{Upper bound of confidence interval}
#'     \item{alpha}{Significance level used}
#'   }
#'
#' @keywords internal
rsv_bootstrap <- function(Y, D, S, pred_Y, pred_D, pred_S_e, pred_S_o,
                          clusters = NULL, B = 1000, alpha = 0.05, eps = 1e-6) {

  # If no clusters provided, use individual-level bootstrap
  if (is.null(clusters)) {
    clusters <- seq_along(Y)
  }

  # Combine data for bootstrap
  boot_data <- data.frame(
    Y = Y,
    D = D,
    S = S,
    pred_Y = pred_Y,
    pred_D = pred_D,
    pred_S_e = pred_S_e,
    pred_S_o = pred_S_o,
    clusters = clusters
  )

  # Bootstrap function
  boot_stat <- function(data, indices) {
    d <- data[indices, ]
    res <- rsv_compute(
      Y = d$Y,
      D = d$D,
      S = d$S,
      pred_Y = d$pred_Y,
      pred_D = d$pred_D,
      pred_S_e = d$pred_S_e,
      pred_S_o = d$pred_S_o,
      eps = eps
    )
    res$coef
  }

  # Cluster bootstrap
  unique_clusters <- unique(clusters)
  boot_coefs <- replicate(B, {
    sampled_clusters <- sample(unique_clusters, length(unique_clusters), replace = TRUE)
    indices <- which(clusters %in% sampled_clusters)
    boot_stat(boot_data, indices)
  })

  # Compute SE and confidence interval
  se <- sd(boot_coefs, na.rm = TRUE)
  ci_lower <- quantile(boot_coefs, alpha / 2, na.rm = TRUE)
  ci_upper <- quantile(boot_coefs, 1 - alpha / 2, na.rm = TRUE)

  list(
    se = se,
    ci_lower = unname(ci_lower),
    ci_upper = unname(ci_upper),
    alpha = alpha
  )
}
