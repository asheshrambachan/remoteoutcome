# =============================================================================
# RSV with Sample Splitting
# =============================================================================

#' RSV Treatment Effect Estimator with Sample Splitting
#'
#' Estimates treatment effects using remotely sensed variables (RSVs) with
#' sample splitting. Splits data into training and test sets, fits prediction
#' models on training set, and evaluates on test set.
#'
#' @inheritParams rsv
#' @param y_levels Levels of Y outcome. Defaults to all unique non-NA values of Y, sorted.
#' @param clusters Optional cluster identifiers. If provided, the train/test split
#'   is performed at the cluster level.
#' @param theta_init Optional initial treatment effect estimate. If \code{NULL},
#'   estimated from training data.
#' @param train.proportion Proportion of data to use for training (default 0.5).
#'   Must be between 0 and 1 (exclusive).
#'
#' @return An object of class \code{"rsv"}. See \code{\link{rsv}} for details.
#'
#' @details
#' This function performs sample splitting to reduce overfitting bias:
#' \enumerate{
#'   \item Randomly splits data into training and test sets
#'   \item Fits prediction models on training set
#'   \item Makes predictions on test set
#'   \item Computes RSV estimate using test set
#' }
#'
#' Sample splitting provides less biased estimates than \code{\link{rsv}} but
#' uses only a subset of data for estimation. For a more efficient approach,
#' consider \code{\link{cv.rsv}} which uses cross-fitting.
#'
#' Standard errors can be added using \code{\link{add_se}}.
#'
#' @seealso \code{\link{rsv}}, \code{\link{cv.rsv}}, \code{\link{rsv_fitted}},
#'   \code{\link{add_se}}
#'
#' @examples
#' \dontrun{
#' # 50/50 train/test split
#' fit <- rsv_split(Y, D, S_e, S_o, R,
#'                  models = list(
#'                    Y = list(model = "rf", num.trees = 500, class.weights = c(10, 1)),
#'                    D = list(model = "rf", num.trees = 500),
#'                    S_e = list(model = "rf", num.trees = 500),
#'                    S_o = list(model = "rf", num.trees = 500)
#'                  ),
#'                  train.proportion = 0.5,
#'                  seed = 42)
#'
#' # 70/30 train/test split with mixed models
#' fit <- rsv_split(Y, D, S_e, S_o, R,
#'                  models = list(
#'                    Y = list(model = "rf", num.trees = 200, class.weights = c(10, 1)),
#'                    D = list(model = "rf", num.trees = 200),
#'                    S_e = list(model = "logit"),
#'                    S_o = list(model = "logit")
#'                  ),
#'                  train.proportion = 0.7,
#'                  seed = 42)
#'
#' # Add standard errors
#' fit <- add_se(fit, B = 1000, clusters = village_id, num.threads = 8)
#' }
#'
#' @export
rsv_split <- function(Y, D, S_e, S_o, R,
                      X = NULL,
                      y_levels = NULL,
                      yK = NULL,
                      models,
                      train.proportion = 0.5,
                      clusters = NULL,
                      sigma2.lower = NULL,
                      sigma2.quantile = 0.01,
                      theta_init = NULL,
                      seed = NULL,
                      num.threads = 1) {

  # Store call
  cl <- match.call()

  # Validate inputs
  validate_required_inputs(Y, D, S_e, S_o)
  n <- length(Y)
  R <- validate_R(R, n)
  validate_split_params(train.proportion)

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

  # Create train/test split (cluster-level if clusters provided)
  split <- make_split(n, train.proportion, seed, clusters)
  train_idx <- split$train
  test_idx <- split$test

  # Fit models on training data, predict on test data
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

  # Get theta_init
  theta_init_val <- theta_init %||% fit_result$theta_init

  # Combine observations, predictions into single df
  df <- cbind(
    data.frame(
      Y = Y[test_idx],
      D = D[test_idx],
      S_e = S_e[test_idx],
      S_o = S_o[test_idx]
    ),
    fit_result$predictions
  )

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
    models = model_specs,
    split = "sample",
    train.proportion = train.proportion,
    test.size = length(test_idx),
    clusters = clusters,
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
