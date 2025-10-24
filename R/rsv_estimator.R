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
#'   \item Provide fitted predictions \code{pred_Y}, \code{pred_D}, \code{pred_S_e}, \code{pred_S_o} directly.
#'   \item Provide raw data \code{(Y, D, S, R)} and the function fits predictions using random forests.
#' }
#'
#' The function also supports sample splitting and K-fold cross-fitting for
#' prediction fitting.
#'
#' @param Y Outcome variable (binary, \code{NA} where not observed). 
#' @param D Treatment indicator (binary, \code{NA} where not observed). 
#' @param S_e Experimental sample indicator (0 or 1).
#' @param S_o Observational sample indicator (0 or 1).
#' @param R Remotely sensed variable. Required if predictions are not provided.
#' @param pred_Y (Optional) Predicted \eqn{E[Y \mid R, S_o = 1]}. If provided, other predictions must also be provided.
#' @param pred_D (Optional) Predicted \eqn{E[D \mid R, S_e = 1]}.
#' @param pred_S_e (Optional) Predicted \eqn{P(S_e = 1 \mid R)}.
#' @param pred_S_o (Optional) Predicted \eqn{P(S_o = 1 \mid R)}.
#' @param theta_init Initial estimate of the treatment effect on the train data.
#' @param eps Small constant for numerical stability of \code{sigma2} estimate (default \code{1e-2}).
#' @param method Prediction fitting method; one of \code{"split"} (default), \code{"crossfit"}, or \code{"none"}.
#'   \code{"split"} = simple sample split; \code{"crossfit"} = K-fold cross-fitting; \code{"none"} = use all data for training/testing.
#' @param ml_params List of parameters for random forest:
#'   \describe{
#'     \item{ntree}{Number of trees}
#'     \item{classwt_Y}{Class weights for \code{pred_Y} model; default \code{c(10, 1)})}
#'     \item{seed}{User specified seed passed to each \code{ranger} function for reproducibility; default \code{NULL}}
#'     \item{cores}{Number of cores used by `ranger` for parallel training; default \code{1}}
#'     \item{nfolds}{Number of folds for cross-fitting (default \code{5}).}
#'     \item{train_ratio}{Proportion for training in sample split (default \code{0.5}).}
#'   }
#' @param se Logical; compute standard errors via bootstrap? (default \code{TRUE}).
#' @param se_params List of bootstrap parameters:
#'   \describe{
#'     \item{B}{Number of bootstrap replications (default \code{1000})}
#'     \item{fix_seed}{If \code{TRUE}, deterministic seeding is used with `set.seed(b)` for the *b*‑th replication (default \code{FALSE})}
#'     \item{cores}{Number of cores for bootstrap replications (default \code{1}).}
#'     \item{clusters}{Clusters for the bootstrap. If \code{NULL}, uses individual-level bootstrap}
#'   }
#'
#' @return A list of class \code{"rsv"} with components:
#' \describe{
#'   \item{coef}{Treatment effect estimate.}
#'   \item{se}{Standard error (if \code{se = TRUE}).}
#'   \item{denominator_se}{Standard error of the denominator of the treatment effect (if \code{se = TRUE})}
#'   \item{n_obs}{Sample size in observational sample.}
#'   \item{n_exp}{Sample size in experimental sample.}
#'   \item{n_both}{Sample size in both samples.}
#'   \item{numerator}{Numerator of the treatment effect estimate}
#'   \item{denominator}{Denominator of the treatment effect estimate}
#'   \item{method}{Prediction fitting method used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: User provides raw data with sample splitting
#' result <- rsv_estimate(Y = Y, D = D, S_e = S_e, S_o = S_o, R = R,
#'                        method = "split", se = TRUE)
#'
#' # Example 2: User provides raw data with cross-fitting
#' result <- rsv_estimate(Y = Y, D = D, S_e = S_e, S_o = S_o, R = R,
#'                        method = "crossfit", nfolds = 5, se = TRUE)
#'
#' # Example 3: User provides raw data with custom ML parameters
#' result <- rsv_estimate(Y = Y, D = D, S_e = S_e, S_o = S_o, R = R,
#'                        method = "crossfit",
#'                        ml_params = list(ntree = 1000, classwt_Y = c(5, 1)),
#'                        se = TRUE)
#'
#' # Example 4: User provides fitted predictions
#' result <- rsv_estimate(Y = Y, D = D, S_e = S_e, S_o = S_o,
#'                        pred_Y = pred_Y, pred_D = pred_D,
#'                        pred_S_e = pred_S_e, pred_S_o = pred_S_o,
#'                        se = TRUE)
#' }
#' 
rsv_estimate <- function(
  Y = NULL, D = NULL, S_e = NULL, S_o = NULL, R = NULL,
  pred_Y = NULL, pred_D = NULL, pred_S_e = NULL, pred_S_o = NULL, theta_init = NULL,
  eps = 1e-2,
  method = c("split", "crossfit", "none"), ml_params = list(),
  se = TRUE, se_params = list()
  ) {
  
  # Match method argument
  method <- match.arg(method)

  # Store the call
  call <- match.call()
  
  # Merge user-provided ML params with defaults
  ml_defaults <- list(
    ntree = 100,
    classwt_Y = c(10, 1),
    seed = NULL,
    cores = 1,
    train_ratio = 0.5, # method = "split" parameter
    nfolds = 5 # method = "crossfit" parameter
  )
  ml_params <- utils::modifyList(ml_defaults, ml_params)

  # Merge user-provided se boot params with defaults
  se_defaults <- list(
    B = 1000,
    fix_seed = FALSE,
    cores = 1,
    clusters = NULL
  )
  se_params <- utils::modifyList(se_defaults, se_params)
  
  # Check if Y, D, S_e, and S_o are provided
  if (is.null(Y) || is.null(D) || is.null(S_e) || is.null(S_o))
    stop("Y, D, S_e, and S_o must be provided")
  
  # Check if predictions are provided
  preds_provided <- !is.null(pred_Y) && !is.null(pred_D) &&
                    !is.null(pred_S_e) && !is.null(pred_S_o)
  
  if (preds_provided) {
    # Interface 1: User provides predictions
    result <- rsv_from_predictions(
      Y, D, S_e, S_o, 
      pred_Y, pred_D, pred_S_e, pred_S_o, theta_init = theta_init,
      eps = eps,
      ml_params = ml_params,
      se = se, se_params = se_params
    )
    method <- "none"
  } else {
    # Interface 2: Fit predictions from raw data
    if (is.null(R))
      stop("R must be provided to fit predictions from raw data")
    
    # Fit predictions based on method
    if (method == "none")
      result <- rsv_fit_none(
        R = R, Y = Y, D = D, S_e = S_e, S_o = S_o,
        eps = eps,
        ml_params = ml_params,
        se = se, se_params = se_params
      )
    else if (method == "split")
      result <- rsv_fit_split(
        R = R, Y = Y, D = D, S_e = S_e, S_o = S_o,
        eps = eps,
        ml_params = ml_params,
        se = se, se_params = se_params
      )
    else if (method == "crossfit")
      result <- rsv_fit_crossfit(
        R = R, Y = Y, D = D, S_e = S_e, S_o = S_o,
        eps = eps,
        ml_params = ml_params,
        se = se, se_params = se_params
      )
  }
  
  # Add metadata
  result$method <- method
  result$call <- call

  # Set class
  class(result) <- "rsv"

  return(result)
}

  
#' Estimate from predictions without fitting
#'
#' Uses all data as a test set
#'
#' @inheritParams rsv_estimate
#' @keywords internal
rsv_from_predictions <- function(
  Y, D, S_e, S_o, 
  pred_Y, pred_D, pred_S_e, pred_S_o, 
  theta_init = NULL,
  eps = 1e-2,
  ml_params = list(),
  se = TRUE,
  se_params = list()
) {
  
  # Construct test set
  if (is.null(theta_init))
    # Estimate theta_init from test data instead of train data if not provided
    theta_init <- get_theta_init(
      observations = data.frame(Y, D, S_e, S_o),
      predictions = data.frame(
        Y = pred_Y,
        D = pred_D,
        S_e = S_e, # pred_S_e,
        S_o = S_o # pred_S_o
      )
    )
  
  predictions <- data.frame(
    Y = pred_Y,
    D = pred_D,
    S_e = pred_S_e,
    S_o = pred_S_o
  )
  observations <- data.frame(Y, D, S_e, S_o)
  
  # Compute RSV estimator
  result <- rsv_compute(
    observations = observations,
    predictions = predictions,
    theta_init = theta_init, 
    eps = eps
  )
  
  # Compute standard errors if requested
  if (se) {
    boot_result <- rsv_bootstrap(
      observations = observations,
      predictions = predictions,
      theta_init = theta_init,
      eps = eps,
      se_params = se_params
    )
    result <- c(result, boot_result)
  }
  
  return(result)
}

#' Fit predictions with no split
#'
#' Uses all data for both training and testing.
#'
#' @inheritParams rsv_estimate
#' @keywords internal
rsv_fit_none <- function(
    R, Y, D, S_e, S_o, 
    eps = 1e-2,
    ml_params = list(),
    se = TRUE,
    se_params = list()
  ) {
  
  # Fit on training data, predict on test data
  out <- fit_predictions_rf(
    R = R, Y = Y, D = D, S_e = S_e, S_o = S_o,
    R_pred = R,
    ml_params = ml_params
  )
  
  # Construct test set
  theta_init <- out$theta_init
  predictions <- out$predictions
  observations <- data.frame(Y, D, S_e, S_o)
  
  # Compute RSV estimator
  result <- rsv_compute(
    observations = observations,
    predictions = predictions,
    theta_init = theta_init, 
    eps = eps
  )
  
  # Compute standard errors if requested
  if (se) {
    boot_result <- rsv_bootstrap(
      observations = observations,
      predictions = predictions,
      theta_init = theta_init,
      eps = eps,
      se_params = se_params
    )
    result <- c(result, boot_result)
  }
  
  return(result)
}

#' Fit predictions with sample splitting
#'
#' Splits data into training and test sets, fits predictions on training set,
#' and returns predictions on test set.
#'
#' @inheritParams rsv_estimate
#' @keywords internal
rsv_fit_split <- function(
    R, Y, D, S_e, S_o, 
    eps = 1e-2,
    ml_params = list(),
    se = TRUE,
    se_params = list()
  ) {
  
  n <- length(Y)
  train_ratio <- ml_params$train_ratio
  seed <- ml_params$seed
  
  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Create train/test split
  train_idx <- sample(seq_len(n), size = floor(train_ratio * n), replace = FALSE)
  test_idx <- setdiff(seq_len(n), train_idx)

  # Fit on training data, predict on test data
  out <- fit_predictions_rf(
    R = R[train_idx, , drop=FALSE],
    Y = Y[train_idx],
    D = D[train_idx],
    S_e = S_e[train_idx],
    S_o = S_o[train_idx],
    R_pred = R[test_idx, , drop=FALSE],
    ml_params = ml_params
  )
  
  # Construct test set
  theta_init <- out$theta_init
  predictions <- out$predictions
  observations <- data.frame(
    Y = Y[test_idx], 
    D = D[test_idx], 
    S_e = S_e[test_idx],
    S_o = S_o[test_idx]
  )
  
  # Compute RSV estimator
  result <- rsv_compute(
    observations = observations,
    predictions = predictions,
    theta_init = theta_init, 
    eps = eps
  )
  
  # Compute standard errors if requested
  if (se) {
    if (!is.null(se_params$clusters))
      se_params$clusters <- se_params$clusters[test_idx]
    
    boot_result <- rsv_bootstrap(
      observations = observations,
      predictions = predictions,
      theta_init = theta_init,
      eps = eps,
      se_params = se_params
    )
    result <- c(result, boot_result)
  }
  
  return(result)
}


#' Fit predictions with K-fold cross-fitting
#'
#' Performs K-fold cross-fitting: splits data into K folds, fits predictions
#' on K-1 folds, predicts on held-out fold, and repeats for all folds.
#'
#' @inheritParams rsv_estimate
#' @keywords internal
rsv_fit_crossfit <- function(
    R, Y, D, S_e, S_o, 
    eps = 1e-2,
    ml_params = list(),
    se = FALSE,
    se_params = list()
  ) {

  n <- length(Y)
  nfolds <- ml_params$nfolds
  seed <- ml_params$seed
  
  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Create folds
  fold_ids <- sample(rep(seq_len(nfolds), length.out = n))

  # Initialize empty lists to store results across folds
  results <- vector("list", nfolds)
  coefs <- rep(NA, nfolds)
  
  # Cross-fitting loop
  for (k in seq_len(nfolds)) {
    # Training set: all folds except k
    train_idx <- (fold_ids != k)
    # Test set: fold k
    test_idx <- (fold_ids == k)

    # Fit on training folds, predict on test fold
    out_k <- fit_predictions_rf(
      R = R[train_idx, , drop = FALSE],
      Y = Y[train_idx],
      D = D[train_idx],
      S_e = S_e[train_idx],
      S_o = S_o[train_idx],
      R_pred = R[test_idx, , drop = FALSE],
      ml_params = ml_params
    )
    
    # Construct test set
    theta_init_k <- out_k$theta_init
    predictions_k <- out_k$predictions
    observations_k <- data.frame(
      Y = Y[test_idx], 
      D = D[test_idx], 
      S_e = S_e[test_idx],
      S_o = S_o[test_idx]
    )
    
    # Compute RSV estimator
    result_k <- rsv_compute(
      observations = observations_k,
      predictions = predictions_k,
      theta_init = theta_init_k, 
      eps = eps
    )
    
    # Compute standard errors if requested
    if (se) {
      se_params_k <- se_params
      
      if (!is.null(se_params_k$clusters))
        se_params_k$clusters <- se_params_k$clusters[test_idx]
      
      boot_result_k <- rsv_bootstrap(
        observations = observations_k,
        predictions = predictions_k,
        theta_init = theta_init_k,
        eps = eps,
        se_params = se_params_k
      )
      result_k <- c(result_k, boot_result_k)
    }
    
    results[[k]] <- result_k
    coefs[k] <- result_k$coef
  }
  
  # Pooled Cross-Validated Estimator
  coef <- mean(coefs, na.rm = TRUE) 
  
  return(list(coef = coef, results = results))
}
