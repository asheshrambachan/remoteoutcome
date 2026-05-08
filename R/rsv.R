# =============================================================================
# RSV Main Function (No Sample Splitting)
# =============================================================================

#' RSV Treatment Effect Estimator (No Sample Splitting)
#'
#' Estimates treatment effects using remotely sensed variables (RSVs) without
#' sample splitting. Fits prediction models on the full dataset and evaluates
#' on the same data (resubstitution).
#'
#' @param Y Outcome variable (can be multi-category, \code{NA} where not observed).
#' @param D Treatment indicator (binary, \code{NA} where not observed).
#' @param S_e Experimental sample indicator (logical or 0/1).
#' @param S_o Observational sample indicator (logical or 0/1).
#' @param R Remote sensing variables (matrix or data frame).
#' @param X Discrete covariates (matrix or data frame, optional). All columns will be used.
#' @param y_levels Levels of Y outcome. Defaults to all unique non-NA values of Y, sorted.
#' @param yK Reference category for Y.
#' @param models Named list providing model specifications for all four models.
#'   Must contain elements \code{Y}, \code{D}, \code{S_e}, \code{S_o}, where each
#'   element is a list specifying the model type and parameters. See Details.
#' @param sigma2.lower Lower bound for sigma-squared stabilization. If \code{NULL}
#'   (default), uses quantile-based lower bound.
#' @param sigma2.quantile Quantile for sigma-squared lower bound (default 0.01).
#'   Only used if \code{sigma2.lower = NULL}. Sigma-squared values are bounded
#'   below by their \code{sigma2.quantile}-th quantile for numerical stability.
#' @param seed Random seed for reproducibility (default \code{NULL}).
#' @param num.threads Number of threads for parallel processing in random forest
#'   (default 1).
#'
#' @return An object of class \code{"rsv"} containing:
#' \describe{
#'   \item{coefficients}{Treatment effect estimate (named numeric)}
#'   \item{numerator}{Numerator of the treatment effect estimate}
#'   \item{denominator}{Denominator of the treatment effect estimate}
#'   \item{theta_init}{Initial treatment effect estimate}
#'   \item{weights}{Efficient weights used in estimation}
#'   \item{n.obs}{Sample size in observational sample}
#'   \item{n.exp}{Sample size in experimental sample}
#'   \item{n.both}{Sample size in both samples}
#'   \item{model}{Model type(s) used}
#'   \item{split}{Character: "none" for this function}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' This function fits prediction models for \eqn{P(Y|R, S_o=1)}, \eqn{P(D|R, S_e=1)},
#' \eqn{P(S_e=1|R)}, and \eqn{P(S_o=1|R)} using the full dataset, then uses
#' those predictions to estimate the treatment effect on the same data. This is
#' known as resubstitution and may lead to overfitting. For less biased estimates,
#' use \code{\link{rsv_split}} or \code{\link{cv.rsv}}.
#'
#' The \code{models} argument allows fine-grained control over each prediction model:
#' \preformatted{
#' models = list(
#'   Y = list(model = "rf", num.trees = 1000, class.weights = c(10, 1)),
#'   D = list(model = "rf", num.trees = 500),
#'   S_e = list(model = "logit"),
#'   S_o = list(model = "logit")
#' )
#' }
#'
#' Standard errors can be added using \code{\link{add_se}}.
#'
#' @seealso \code{\link{rsv_split}}, \code{\link{cv.rsv}}, \code{\link{rsv_fitted}},
#'   \code{\link{add_se}}
#'
#' @examples
#' \dontrun{
#' # All models using random forest
#' fit <- rsv(Y, D, S_e, S_o, R,
#'            models = list(
#'              Y = list(model = "rf", num.trees = 500, class.weights = c(10, 1)),
#'              D = list(model = "rf", num.trees = 500),
#'              S_e = list(model = "rf", num.trees = 500),
#'              S_o = list(model = "rf", num.trees = 500)
#'            ),
#'            seed = 42)
#'
#' # Mixed models (RF for Y and D, logistic for S_e and S_o)
#' fit <- rsv(Y, D, S_e, S_o, R,
#'            models = list(
#'              Y = list(model = "rf", num.trees = 1000, class.weights = c(10, 1)),
#'              D = list(model = "rf", num.trees = 500),
#'              S_e = list(model = "logit"),
#'              S_o = list(model = "logit")
#'            ),
#'            seed = 42)
#'
#' # Add standard errors
#' fit <- add_se(fit, B = 1000, clusters = village_id, num.threads = 8)
#' }
#'
#' @export
rsv <- function(Y, D, S_e, S_o, R,
                X = NULL,
                y_levels,
                yK = NULL,
                models,
                sigma2.lower = NULL,
                sigma2.quantile = 0.01,
                seed = NULL,
                num.threads = 1) {

  # Store call
  cl <- match.call()

  # Validate required inputs
  validate_required_inputs(Y, D, S_e, S_o)
  n <- length(Y)
  R <- validate_R(R, n)

  # Validate X if provided and infer X_cols
  X_cols <- NULL
  if (!is.null(X)) {
    X <- validate_X(X, n, NULL)
    X_cols <- colnames(X)
  }

  # Resolve model specifications
  model_specs <- resolve_models(models)
  validate_models(model_specs)

  # Fit models on full data, predict on full data (resubstitution)
  fit_result <- fit_all_models(
    R_train = R,
    X_train = X,
    Y = Y,
    D = D,
    S_e = S_e,
    S_o = S_o,
    R_test = R,
    X_test = X,
    X_cols = X_cols,
    y_levels = y_levels,
    yK = yK,
    models = model_specs,
    seed = seed,
    num.threads = num.threads
  )

  # Get theta_init
  theta_init_val <- fit_result$theta_init

  # Combine observations, predictions into single df
  df <- cbind(
    data.frame(Y = Y, D = D, S_e = S_e, S_o = S_o),
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
    split = "none",
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
