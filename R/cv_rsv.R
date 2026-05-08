# =============================================================================
# RSV with K-Fold Cross-Fitting
# =============================================================================
#' RSV Treatment Effect Estimator with Cross-Fitting
#'
#' Estimates treatment effects using remotely sensed variables (RSVs) with
#' K-fold cross-fitting. Provides more efficient use of data than simple
#' sample splitting.
#'
#' @inheritParams rsv
#' @param y_levels Levels of Y outcome (default: all unique non-NA values of Y, sorted).
#' @param clusters Optional cluster identifiers for clustered cross-fitting (default \code{NULL}).
#' @param nfolds Number of folds for cross-validation (default 5). Must be
#'   at least 2 and cannot exceed sample size.
#'
#' @return An object of class \code{"cv.rsv"} containing:
#' \describe{
#'   \item{coefficients}{Pooled treatment effect estimate (single compute_estimate call on all n obs with cross-fitted predictions)}
#'   \item{relevance}{ATE-specific relevance \eqn{(\lambda' J^{-1} \lambda)^{-1}}}
#'   \item{coefficients_naive}{Naive ATE estimate (H = P(Y=yj|D=d,R))}
#'   \item{relevance_naive}{Naive estimator relevance}
#'   \item{weights}{Observation weights from the pooled estimate}
#'   \item{nfolds}{Number of folds used}
#'   \item{n.obs}{Total sample size in observational sample}
#'   \item{n.exp}{Total sample size in experimental sample}
#'   \item{n.both}{Total sample size in both samples}
#'   \item{model}{Model type(s) used}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' This function performs K-fold cross-fitting:
#' \enumerate{
#'   \item Randomly splits data into K folds
#'   \item For each fold k:
#'     \itemize{
#'       \item Fits prediction models on all folds except k
#'       \item Makes predictions on fold k
#'       \item Computes RSV estimate on fold k
#'     }
#'   \item Returns average estimate across all folds
#' }
#'
#' Cross-fitting provides less biased estimates than resubstitution
#' (\code{\link{rsv}}) and uses data more efficiently than simple sample
#' splitting (\code{\link{rsv_split}}).
#'
#' Standard errors can be added using \code{\link{add_se}}, which will
#' bootstrap the entire cross-fitting procedure.
#'
#' @seealso \code{\link{rsv}}, \code{\link{rsv_split}}, \code{\link{rsv_fitted}},
#'   \code{\link{add_se}}
#'
#' @examples
#' \dontrun{
#' # 5-fold cross-fitting with random forest
#' fit <- cv.rsv(Y, D, S_e, S_o, R,
#'               models = list(
#'                 Y = list(model = "rf", num.trees = 500, class.weights = c(10, 1)),
#'                 D = list(model = "rf", num.trees = 500),
#'                 S_e = list(model = "rf", num.trees = 500),
#'                 S_o = list(model = "rf", num.trees = 500)
#'               ),
#'               nfolds = 5,
#'               seed = 42)
#'
#' # 10-fold cross-fitting with mixed models
#' fit <- cv.rsv(Y, D, S_e, S_o, R,
#'               models = list(
#'                 Y = list(model = "rf", num.trees = 1000, class.weights = c(10, 1)),
#'                 D = list(model = "rf", num.trees = 500),
#'                 S_e = list(model = "logit"),
#'                 S_o = list(model = "logit")
#'               ),
#'               nfolds = 10,
#'               seed = 42)
#'
#' # Examine fold-specific estimates
#' print(fit)
#' summary(fit)
#' plot(fit)
#'
#' # Add standard errors (expensive - bootstraps entire CV)
#' fit <- add_se(fit, B = 500, clusters = village_id, num.threads = 8)
#' }
#'
#' @export
cv.rsv <- function(Y, D, S_e, S_o, R,
                   X = NULL,
                   y_levels = NULL,
                   yK = NULL,
                   models,
                   nfolds = 5,
                   clusters = NULL,
                   sigma2.lower = NULL,
                   sigma2.quantile = 0.01,
                   seed = NULL,
                   num.threads = 1) {

  # Store call
  cl <- match.call()

  # Validate inputs
  validate_required_inputs(Y, D, S_e, S_o)
  n <- length(Y)
  R <- validate_R(R, n)
  validate_cv_params(nfolds, n)

  if (is.null(y_levels))
    y_levels <- sort(unique(na.omit(Y)))
  
  # Validate X if provided and infer X_cols
  X_cols <- NULL
  if (!is.null(X)) {
    X <- validate_X(X, n, NULL)
    X_cols <- colnames(X)
  }

  # Resolve model specifications
  model_specs <- resolve_models(models)
  validate_models(model_specs)

  # Create folds (cluster-level if clusters provided)
  fold_ids <- make_folds(n, nfolds, seed, clusters)
  
  # Storage for pooled cross-fitted predictions
  pooled_preds_list   <- vector("list", nfolds)
  pooled_test_indices <- vector("list", nfolds)

  # Loop over folds: fit models on training folds, predict on test fold
  for (k in seq_len(nfolds)) {
    train_idx <- (fold_ids != k)
    test_idx  <- (fold_ids == k)

    fit_result <- fit_all_models(
      R_train = R[train_idx, , drop = FALSE],
      X_train = if (!is.null(X)) X[train_idx, , drop = FALSE] else NULL,
      Y = Y[train_idx],
      D = D[train_idx],
      S_e = S_e[train_idx],
      S_o = S_o[train_idx],
      R_test = R[test_idx, , drop = FALSE],
      X_test = if (!is.null(X)) X[test_idx, , drop = FALSE] else NULL,
      X_cols = X_cols,
      y_levels = y_levels,
      yK = yK,
      models = model_specs,
      seed = seed,
      num.threads = num.threads
    )

    pooled_preds_list[[k]]   <- fit_result$predictions
    pooled_test_indices[[k]] <- which(test_idx)
  }

  # Assemble pooled cross-fitted predictions in original observation order
  all_test_idx <- unlist(pooled_test_indices)
  all_preds <- do.call(rbind, pooled_preds_list)
  pooled_predictions <- all_preds[order(all_test_idx), , drop = FALSE]
  rownames(pooled_predictions) <- NULL

  # Pooled estimate: single compute_estimate call on all n obs with cross-fitted
  # predictions. This is the estimator the tex theory targets and that the SE
  # methods bootstrap — solving (1/n) sum_i G_i(theta, h_{k(i)}) = 0.
  yK_use    <- if (is.null(yK)) y_levels[1] else yK
  X_cols_use <- if (!is.null(X_cols) && length(X_cols) > 0) X_cols else character(0)
  df_pooled <- cbind(
    data.frame(Y = Y, D = D, S_e = S_e, S_o = S_o),
    pooled_predictions
  )
  df_pooled <- .ensure_numeric_indicators(df_pooled)
  theta_init_pooled <- compute_theta_init(df_pooled, y_levels = y_levels,
                                          yK = yK_use, X_cols = X_cols_use)
  pooled_result <- compute_estimate(
    df              = df_pooled,
    theta_init      = theta_init_pooled,
    y_levels        = y_levels,
    sigma2.lower    = sigma2.lower,
    sigma2.quantile = sigma2.quantile
  )

  relevance <- pooled_result$relevance_ate

  # Naive estimate: H = P(Y=yj|D=d,R) (always computed alongside RSV)
  naive_result <- compute_estimate(
    df              = df_pooled,
    theta_init      = theta_init_pooled,
    y_levels        = y_levels,
    sigma2.lower    = sigma2.lower,
    sigma2.quantile = sigma2.quantile,
    pred_y_only     = TRUE
  )
  relevance_naive <- naive_result$relevance_ate

  # Total sample sizes
  n.obs  <- sum(S_o == 1 | S_o == TRUE)
  n.exp  <- sum(S_e == 1 | S_e == TRUE)
  n.both <- sum((S_e == 1 | S_e == TRUE) & (S_o == 1 | S_o == TRUE))

  # Construct result object
  result <- list(
    coefficients       = pooled_result$coefficients,
    relevance          = relevance,
    coefficients_naive = naive_result$coefficients,
    relevance_naive    = relevance_naive,
    weights            = pooled_result$weights,
    weights_naive      = naive_result$weights,
    nfolds = nfolds,
    n.obs = n.obs,
    n.exp = n.exp,
    n.both = n.both,
    models = model_specs,
    call = cl,
    # Store for bootstrap (entire dataset needed for CV bootstrap)
    Y = Y,
    D = D,
    S_e = S_e,
    S_o = S_o,
    R = R,
    X = X,
    X_cols = X_cols,
    y_levels = y_levels,
    yK = yK,
    clusters = clusters,
    fold_ids = fold_ids,
    pooled_predictions = pooled_predictions,
    sigma2.lower = sigma2.lower,
    sigma2.quantile = sigma2.quantile
  )

  class(result) <- "cv.rsv"
  return(result)
}
