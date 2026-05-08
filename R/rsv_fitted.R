# =============================================================================
# RSV from Pre-fitted Predictions
# =============================================================================

#' RSV Treatment Effect Estimator from Pre-fitted Predictions
#'
#' Estimates treatment effects using user-provided predicted probabilities.
#' Useful when you have already fitted prediction models externally or want
#' to use custom prediction methods.
#'
#' @param Y Outcome variable (can be multi-category, \code{NA} where not observed).
#' @param D Treatment indicator (binary, \code{NA} where not observed).
#' @param S_e Experimental sample indicator (logical or 0/1).
#' @param S_o Observational sample indicator (logical or 0/1).
#' @param X Discrete covariates (matrix or data frame, optional). All columns will be used.
#' @param y_levels Levels of Y outcome (default c(0, 1) for binary).
#' @param yK Reference category for Y (default first level in y_levels).
#' @param pred_Y Matrix or data.frame with predicted probabilities for all Y levels.
#'   Columns should be named pred_Y0, pred_Y1, pred_Y2, ... matching y_levels.
#'   For binary Y, can also provide a vector of P(Y=1).
#' @param pred_D Predicted \eqn{P(D=1 | X, R, S_e=1)}.
#' @param pred_S_e Predicted \eqn{P(S_e=1 | X, R)}.
#' @param pred_S_o Predicted \eqn{P(S_o=1 | X, R)}.
#' @param theta_init Initial estimate of treatment effect (optional). If \code{NULL},
#'   will be estimated from a subset of the data.
#' @param train.proportion Proportion of data to use for estimating \code{theta_init}
#'   if it is not provided (default 0.2). Only used when \code{theta_init = NULL}.
#' @param sigma2.lower Lower bound for sigma-squared stabilization. If \code{NULL}
#'   (default), uses quantile-based lower bound.
#' @param sigma2.quantile Quantile for sigma-squared lower bound (default 0.01).
#' @param seed Random seed for reproducibility when estimating \code{theta_init}
#'   (default \code{NULL}).
#'
#' @return An object of class \code{"rsv"}. See \code{\link{rsv}} for details.
#'
#' @details
#' This function allows users to provide their own predicted probabilities,
#' which can come from any source:
#' \itemize{
#'   \item Models fitted using other packages (e.g., xgboost, keras)
#'   \item Ensemble predictions
#'   \item Pre-computed predictions from previous runs
#'   \item Custom machine learning pipelines
#' }
#'
#' If \code{theta_init} is not provided, the function randomly selects
#' \code{train.proportion} of the data to estimate it, then uses the remaining
#' data for the final RSV estimate.
#'
#' Standard errors can be added using \code{\link{add_se}}.
#'
#' @seealso \code{\link{rsv}}, \code{\link{rsv_split}}, \code{\link{cv.rsv}},
#'   \code{\link{add_se}}
#'
#' @examples
#' \dontrun{
#' # With theta_init provided
#' fit <- rsv_fitted(Y, D, S_e, S_o,
#'                   pred_Y = my_pred_Y,
#'                   pred_D = my_pred_D,
#'                   pred_S_e = my_pred_S_e,
#'                   pred_S_o = my_pred_S_o,
#'                   theta_init = 0.05)
#'
#' # Without theta_init (will estimate from 20% of data)
#' fit <- rsv_fitted(Y, D, S_e, S_o,
#'                   pred_Y = my_pred_Y,
#'                   pred_D = my_pred_D,
#'                   pred_S_e = my_pred_S_e,
#'                   pred_S_o = my_pred_S_o,
#'                   train.proportion = 0.2,
#'                   seed = 42)
#'
#' # Add standard errors
#' fit <- add_se(fit, B = 1000, clusters = village_id, num.threads = 8)
#' }
#'
#' @export
rsv_fitted <- function(Y, D, S_e, S_o,
                       X = NULL,
                       y_levels = c(0, 1),
                       yK = y_levels[1],
                       pred_Y, pred_D, pred_S_e, pred_S_o,
                       theta_init = NULL,
                       train.proportion = 0.2,
                       sigma2.lower = NULL,
                       sigma2.quantile = 0.01,
                       seed = NULL) {

  # Store call
  cl <- match.call()

  # Validate inputs
  validate_required_inputs(Y, D, S_e, S_o)
  n <- length(Y)

  # Validate X if provided and infer X_cols
  X_cols <- NULL
  if (!is.null(X)) {
    X <- validate_X(X, n, NULL)
    X_cols <- colnames(X)
  }

  # Convert pred_Y to data.frame with pred_Y0, pred_Y1, ... if needed
  if (is.vector(pred_Y) || (is.matrix(pred_Y) && ncol(pred_Y) == 1)) {
    # Binary case: pred_Y is P(Y=1)
    if (length(y_levels) != 2) {
      stop("For multi-category Y, pred_Y must be a matrix/data.frame with columns for each level")
    }
    pred_Y_df <- data.frame(
      pred_Y0 = 1 - pred_Y,
      pred_Y1 = pred_Y
    )
    colnames(pred_Y_df) <- paste0("pred_Y", y_levels)
  } else {
    # Multi-category case: check that pred_Y has columns for each level
    pred_Y_df <- as.data.frame(pred_Y)
    expected_cols <- paste0("pred_Y", y_levels)
    if (!all(expected_cols %in% colnames(pred_Y_df))) {
      stop("pred_Y must have columns: ", paste(expected_cols, collapse = ", "))
    }
    pred_Y_df <- pred_Y_df[, expected_cols, drop = FALSE]
  }

  # Validate other predictions
  if (length(pred_D) != n) stop("pred_D must have length n")
  if (length(pred_S_e) != n) stop("pred_S_e must have length n")
  if (length(pred_S_o) != n) stop("pred_S_o must have length n")

  # If theta_init not provided, estimate it from a subset
  if (is.null(theta_init)) {

    validate_split_params(train.proportion)

    # Create split for theta_init estimation
    split <- make_split(n, train.proportion, seed)
    train_idx <- split$train
    test_idx <- split$test

    # Estimate theta_init on train set
    train_df <- cbind(
      data.frame(
        Y = Y[train_idx],
        D = D[train_idx],
        S_e = S_e[train_idx],
        S_o = S_o[train_idx]
      ),
      pred_Y_df[train_idx, , drop = FALSE],
      data.frame(
        pred_D = pred_D[train_idx],
        pred_S_e = pred_S_e[train_idx],
        pred_S_o = pred_S_o[train_idx]
      )
    )

    if (!is.null(X) && !is.null(X_cols) && length(X_cols) > 0) {
      train_df <- cbind(X[train_idx, , drop = FALSE], train_df)
    }

    theta_init_val <- compute_theta_init(train_df, X_cols = X_cols, yK = yK)

    # Use test set for final estimate
    df <- cbind(
      data.frame(
        Y = Y[test_idx],
        D = D[test_idx],
        S_e = S_e[test_idx],
        S_o = S_o[test_idx]
      ),
      pred_Y_df[test_idx, , drop = FALSE],
      data.frame(
        pred_D = pred_D[test_idx],
        pred_S_e = pred_S_e[test_idx],
        pred_S_o = pred_S_o[test_idx]
      )
    )

    if (!is.null(X) && !is.null(X_cols) && length(X_cols) > 0) {
      df <- cbind(X[test_idx, , drop = FALSE], df)
    }

  } else {

    # Use all data
    theta_init_val <- theta_init
    df <- cbind(
      data.frame(Y = Y, D = D, S_e = S_e, S_o = S_o),
      pred_Y_df,
      data.frame(
        pred_D = pred_D,
        pred_S_e = pred_S_e,
        pred_S_o = pred_S_o
      )
    )

    if (!is.null(X) && !is.null(X_cols) && length(X_cols) > 0) {
      df <- cbind(X, df)
    }

    test_idx <- seq_len(n)
  }

  # Compute RSV estimate
  estimate <- compute_estimate(
    df = df,
    theta_init = theta_init_val,
    y_levels = y_levels,
    sigma2.lower = sigma2.lower,
    sigma2.quantile = sigma2.quantile
  )

  # Naive estimate: H = P(Y=yj|D=d,R) (always computed alongside RSV)
  naive_result <- compute_estimate(
    df              = df,
    theta_init      = theta_init_val,
    y_levels        = y_levels,
    sigma2.lower    = sigma2.lower,
    sigma2.quantile = sigma2.quantile,
    pred_y_only     = TRUE
  )
  relevance       <- estimate$relevance_ate
  relevance_naive <- naive_result$relevance_ate

  # Construct result object
  result <- list(
    coefficients       = estimate$coefficients,
    numerator          = estimate$numerator,
    relevance          = relevance,
    theta_init         = estimate$theta_init,
    weights            = estimate$weights,
    coefficients_naive = naive_result$coefficients,
    relevance_naive    = relevance_naive,
    weights_naive      = naive_result$weights,
    n.obs = estimate$n.obs,
    n.exp = estimate$n.exp,
    n.both = estimate$n.both,
    models = "user-provided",
    split = "fitted",
    test.size = length(test_idx),
    call = cl,
    # Store combined df for bootstrap
    data = df,
    X_cols = X_cols,
    y_levels = y_levels,
    yK = yK,
    sigma2.lower = sigma2.lower,
    sigma2.quantile = sigma2.quantile
  )

  class(result) <- "rsv"
  return(result)
}
