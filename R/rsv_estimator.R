# =============================================================================
# RSV Estimator Main Function
# =============================================================================

#' RSV Treatment Effect Estimator
#'
#' Estimates treatment effects using remotely sensed variables (RSVs) following
#' Rambachan, Singh, and Viviano (2025). Implements Algorithm 1 from the main
#' text for binary outcomes without pretreatment covariates.
#'
#' The function supports two interfaces:
#' \enumerate{
#'   \item Provide fitted predictions \code{pred_Y}, \code{pred_D}, \code{pred_S_e},
#'         \code{pred_S_o} directly
#'   \item Provide raw data \code{(Y, D, S, R)} and the function fits predictions
#'         using random forests
#' }
#'
#' The function also supports sample splitting and K-fold cross-fitting for
#' prediction fitting.
#'
#' @param Y Outcome variable (binary, NA where not observed). Required if
#'   predictions are not provided.
#' @param D Treatment indicator (binary, NA where not observed). Required if
#'   predictions are not provided.
#' @param S Sample indicator. Can be:
#'   \describe{
#'     \item{"e"}{Experimental sample only (observe D, R)}
#'     \item{"o"}{Observational sample only (observe Y, R)}
#'     \item{"both"}{Both samples (observe Y, D, R)}
#'   }
#' @param R Remotely sensed variable. Required if predictions are not provided.
#' @param pred_Y (Optional) Predicted E[Y | R, S_o = 1]. If provided, other
#'   predictions must also be provided.
#' @param pred_D (Optional) Predicted E[D | R, S_e = 1].
#' @param pred_S_e (Optional) Predicted P(S_e = 1 | R).
#' @param pred_S_o (Optional) Predicted P(S_o = 1 | R).
#' @param method Prediction fitting method. Options:
#'   \describe{
#'     \item{"split"}{Simple sample split (default)}
#'     \item{"crossfit"}{K-fold cross-fitting}
#'     \item{"none"}{No sample splitting (use all data for training and testing)}
#'   }
#' @param nfolds Number of folds for cross-fitting. Default is 5. Only used if
#'   \code{method = "crossfit"}.
#' @param split_ratio Proportion of data to use for training in sample split.
#'   Default is 0.5. Only used if \code{method = "split"}.
#' @param ml_params List of parameters for random forest. Supported parameters:
#'   \describe{
#'     \item{ntree}{Number of trees. Default is 500.}
#'     \item{classwt_Y}{Class weights for outcome model. Default is c(10, 1).}
#'   }
#'   Additional parameters will be passed to \code{randomForest::randomForest}.
#' @param se Logical. Compute standard errors via bootstrap? Default is TRUE.
#' @param clusters Cluster identifiers for clustered bootstrap. If NULL,
#'   individual-level bootstrap is used.
#' @param B Number of bootstrap replications. Default is 1000.
#' @param alpha Significance level for confidence interval. Default is 0.05
#'   for 95% CI. Use 0.10 for 90% CI, 0.01 for 99% CI, etc.
#' @param eps Small constant for numerical stability. Default is 1e-6.
#' @param seed Random seed for reproducibility in sample splitting/cross-fitting.
#'
#' @return An object of class "rsv" containing:
#'   \describe{
#'     \item{coef}{Treatment effect estimate}
#'     \item{se}{Standard error (if \code{se = TRUE})}
#'     \item{ci_lower}{Lower bound of confidence interval (if \code{se = TRUE})}
#'     \item{ci_upper}{Upper bound of confidence interval (if \code{se = TRUE})}
#'     \item{alpha}{Significance level used for CI (if \code{se = TRUE})}
#'     \item{n_obs}{Sample size in observational sample}
#'     \item{n_exp}{Sample size in experimental sample}
#'     \item{n_both}{Sample size in both samples}
#'     \item{method}{Prediction fitting method used}
#'     \item{call}{The matched call}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: User provides raw data with sample splitting
#' result <- rsv_estimate(Y = Y, D = D, S = S, R = R,
#'                        method = "split", se = TRUE)
#'
#' # Example 2: User provides raw data with cross-fitting
#' result <- rsv_estimate(Y = Y, D = D, S = S, R = R,
#'                        method = "crossfit", nfolds = 5, se = TRUE)
#'
#' # Example 3: User provides raw data with custom ML parameters
#' result <- rsv_estimate(Y = Y, D = D, S = S, R = R,
#'                        method = "crossfit",
#'                        ml_params = list(ntree = 1000, classwt_Y = c(5, 1)),
#'                        se = TRUE)
#'
#' # Example 4: User provides fitted predictions
#' result <- rsv_estimate(Y = Y, D = D, S = S,
#'                        pred_Y = pred_Y, pred_D = pred_D,
#'                        pred_S_e = pred_S_e, pred_S_o = pred_S_o,
#'                        se = TRUE)
#' }
rsv_estimate <- function(Y = NULL, D = NULL, S = NULL, R = NULL,
                         pred_Y = NULL, pred_D = NULL,
                         pred_S_e = NULL, pred_S_o = NULL,
                         method = c("split", "crossfit", "none"),
                         nfolds = 5,
                         split_ratio = 0.5,
                         ml_params = list(),
                         se = TRUE,
                         clusters = NULL,
                         B = 1000,
                         alpha = 0.05,
                         eps = 1e-6,
                         seed = NULL) {

  # Match method argument
  method <- match.arg(method)

  # Store the call
  call <- match.call()

  # Set default ML parameters
  ml_defaults <- list(
    ntree = 500,
    classwt_Y = c(10, 1)
  )

  # Merge user-provided ML params with defaults
  ml_params <- utils::modifyList(ml_defaults, ml_params)

  # Check if predictions are provided
  preds_provided <- !is.null(pred_Y) && !is.null(pred_D) &&
                    !is.null(pred_S_e) && !is.null(pred_S_o)

  # Validate inputs
  if (preds_provided) {
    # Interface 1: User provides predictions
    if (is.null(Y) || is.null(D) || is.null(S)) {
      stop("Y, D, and S must be provided along with predictions")
    }
    if (length(pred_Y) != length(Y)) {
      stop("pred_Y must have same length as Y")
    }
    # Use provided predictions
    predictions <- data.frame(
      pred_Y = pred_Y,
      pred_D = pred_D,
      pred_S_e = pred_S_e,
      pred_S_o = pred_S_o
    )
  } else {
    # Interface 2: Fit predictions from raw data
    if (is.null(Y) || is.null(D) || is.null(S) || is.null(R)) {
      stop("Y, D, S, and R must be provided if predictions are not supplied")
    }

    # Fit predictions based on method
    if (method == "split") {
      predictions <- rsv_fit_split(
        R = R, Y = Y, D = D, S = S,
        split_ratio = split_ratio,
        ml_params = ml_params,
        seed = seed
      )
    } else if (method == "crossfit") {
      predictions <- rsv_fit_crossfit(
        R = R, Y = Y, D = D, S = S,
        nfolds = nfolds,
        ml_params = ml_params,
        seed = seed
      )
    } else if (method == "none") {
      predictions <- fit_predictions_rf(
        R = R, Y = Y, D = D, S = S,
        R_pred = R,
        ml_params = ml_params
      )
    }
  }

  # Compute RSV estimator
  results <- rsv_compute(
    Y = Y,
    D = D,
    S = S,
    pred_Y = predictions$pred_Y,
    pred_D = predictions$pred_D,
    pred_S_e = predictions$pred_S_e,
    pred_S_o = predictions$pred_S_o,
    eps = eps
  )

  # Compute standard errors if requested
  if (se) {
    boot_results <- rsv_bootstrap(
      Y = Y,
      D = D,
      S = S,
      pred_Y = predictions$pred_Y,
      pred_D = predictions$pred_D,
      pred_S_e = predictions$pred_S_e,
      pred_S_o = predictions$pred_S_o,
      clusters = clusters,
      B = B,
      alpha = alpha,
      eps = eps
    )
    results$se <- boot_results$se
    results$ci_lower <- boot_results$ci_lower
    results$ci_upper <- boot_results$ci_upper
    results$alpha <- boot_results$alpha
  }

  # Add metadata
  results$method <- method
  results$call <- call

  # Set class
  class(results) <- "rsv"

  return(results)
}


#' Fit predictions with sample splitting
#'
#' Splits data into training and test sets, fits predictions on training set,
#' and returns predictions on test set.
#'
#' @inheritParams rsv_estimate
#' @keywords internal
rsv_fit_split <- function(R, Y, D, S, split_ratio = 0.5,
                          ml_params = list(), seed = NULL) {

  n <- length(R)

  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Create train/test split
  train_idx <- sample(seq_len(n), size = floor(split_ratio * n), replace = FALSE)
  test_idx <- setdiff(seq_len(n), train_idx)

  # Fit on training data, predict on test data
  predictions_test <- fit_predictions_rf(
    R = R[train_idx],
    Y = Y[train_idx],
    D = D[train_idx],
    S = S[train_idx],
    R_pred = R[test_idx],
    ml_params = ml_params
  )

  # Initialize full predictions with NA
  predictions <- data.frame(
    pred_Y = rep(NA, n),
    pred_D = rep(NA, n),
    pred_S_e = rep(NA, n),
    pred_S_o = rep(NA, n)
  )

  # Fill in test set predictions
  predictions[test_idx, ] <- predictions_test

  return(predictions)
}


#' Fit predictions with K-fold cross-fitting
#'
#' Performs K-fold cross-fitting: splits data into K folds, fits predictions
#' on K-1 folds, predicts on held-out fold, and repeats for all folds.
#'
#' @inheritParams rsv_estimate
#' @keywords internal
rsv_fit_crossfit <- function(R, Y, D, S, nfolds = 5,
                             ml_params = list(), seed = NULL) {

  n <- length(R)

  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Create folds
  fold_ids <- sample(rep(seq_len(nfolds), length.out = n))

  # Initialize predictions
  predictions <- data.frame(
    pred_Y = rep(NA, n),
    pred_D = rep(NA, n),
    pred_S_e = rep(NA, n),
    pred_S_o = rep(NA, n)
  )

  # Cross-fitting loop
  for (k in seq_len(nfolds)) {
    # Training set: all folds except k
    train_idx <- (fold_ids != k)
    # Test set: fold k
    test_idx <- (fold_ids == k)

    # Fit on training folds, predict on test fold
    predictions_k <- fit_predictions_rf(
      R = R[train_idx],
      Y = Y[train_idx],
      D = D[train_idx],
      S = S[train_idx],
      R_pred = R[test_idx],
      ml_params = ml_params
    )

    # Store predictions for fold k
    predictions[test_idx, ] <- predictions_k
  }

  return(predictions)
}
