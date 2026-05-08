# =============================================================================
# Utility Functions
# Validation, parameter processing, and helper functions
# =============================================================================

# Default value operator — NULL-safe; returns b when a is NULL.
#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


#' Validate required inputs
#'
#' @param Y Outcome vector
#' @param D Treatment vector
#' @param S_e Experimental sample indicator
#' @param S_o Observational sample indicator
#'
#' @keywords internal
validate_required_inputs <- function(Y, D, S_e, S_o) {
  if (is.null(Y)) stop("Y must be provided")
  if (is.null(D)) stop("D must be provided")
  if (is.null(S_e)) stop("S_e must be provided")
  if (is.null(S_o)) stop("S_o must be provided")

  n <- length(Y)
  if (length(D) != n) stop("D must have same length as Y")
  if (length(S_e) != n) stop("S_e must have same length as Y")
  if (length(S_o) != n) stop("S_o must have same length as Y")

  invisible(TRUE)
}


#' Validate remote sensing variables
#'
#' @param R Remote sensing variables (matrix or data frame)
#' @param n Expected number of observations
#'
#' @keywords internal
validate_R <- function(R, n) {
  if (is.null(R)) {
    stop("R must be provided")
  }

  R <- as.matrix(R)

  if (nrow(R) != n) {
    stop("R must have same number of rows as Y")
  }

  if (ncol(R) == 0) {
    stop("R must have at least one column")
  }

  # Ensure column names exist so cbind() does not use variable names
  # (e.g. "R_train" vs "R_test"), which breaks nnet::multinom prediction.
  if (is.null(colnames(R))) {
    colnames(R) <- paste0("R", seq_len(ncol(R)))
  }

  invisible(R)
}


#' Validate discrete covariates
#'
#' @param X Discrete covariates (matrix or data frame)
#' @param n Expected number of observations
#' @param X_cols Column names for X (optional)
#'
#' @keywords internal
validate_X <- function(X, n, X_cols = NULL) {
  if (is.null(X)) {
    return(NULL)
  }

  X <- as.data.frame(X)

  if (nrow(X) != n) {
    stop("X must have same number of rows as Y")
  }

  if (ncol(X) == 0) {
    stop("X must have at least one column")
  }

  # If X_cols provided, check that they exist in X
  if (!is.null(X_cols) && length(X_cols) > 0) {
    if (!all(X_cols %in% colnames(X))) {
      missing_cols <- setdiff(X_cols, colnames(X))
      stop("X_cols specifies columns not in X: ", paste(missing_cols, collapse = ", "))
    }
    # Return only the specified columns
    X <- X[, X_cols, drop = FALSE]
  }

  invisible(X)
}


#' Validate predictions
#'
#' @param pred_Y Predicted Y
#' @param pred_D Predicted D
#' @param pred_S_e Predicted S_e
#' @param pred_S_o Predicted S_o
#' @param n Expected length
#'
#' @keywords internal
validate_predictions <- function(pred_Y, pred_D, pred_S_e, pred_S_o, n) {
  if (is.null(pred_Y) || is.null(pred_D) || is.null(pred_S_e) || is.null(pred_S_o)) {
    stop("All predictions (pred_Y, pred_D, pred_S_e, pred_S_o) must be provided")
  }

  if (length(pred_Y) != n) stop("pred_Y must have same length as Y")
  if (length(pred_D) != n) stop("pred_D must have same length as Y")
  if (length(pred_S_e) != n) stop("pred_S_e must have same length as Y")
  if (length(pred_S_o) != n) stop("pred_S_o must have same length as Y")

  invisible(TRUE)
}


#' Resolve model specifications
#'
#' Validates and fills in defaults for model specifications. Each model
#' must specify at minimum the model type ("rf" or "logit").
#'
#' @param models Named list with Y, D, S_e, S_o model specifications.
#'   Each element should be a list with at least "model" specified.
#'
#' @return List with Y, D, S_e, S_o model specifications with defaults filled
#'
#' @keywords internal
resolve_models <- function(models) {

  if (is.null(models)) {
    stop("models must be provided. Specify model type and parameters for Y, D, S_e, S_o.\n",
         "Example: models = list(Y = list(model = 'rf', num.trees = 500))")
  }

  # Default values for RF models
  rf_defaults <- list(
    num.trees = 500,
    mtry = NULL
  )

  model_names <- c("Y", "D", "S_e", "S_o")

  # Check all required models are specified
  missing_models <- setdiff(model_names, names(models))
  if (length(missing_models) > 0) {
    stop("models must specify all four models: Y, D, S_e, S_o.\n",
         "Missing: ", paste(missing_models, collapse = ", "))
  }

  # Fill in defaults for each model
  for (nm in model_names) {

    if (!"model" %in% names(models[[nm]])) {
      stop(sprintf("Model specification for '%s' must include 'model' (either 'rf' or 'logit')", nm))
    }

    # Fill in RF defaults if not specified
    if (models[[nm]]$model == "rf") {
      if (!"num.trees" %in% names(models[[nm]])) {
        models[[nm]]$num.trees <- rf_defaults$num.trees
      }
      if (!"mtry" %in% names(models[[nm]])) {
        models[[nm]]$mtry <- rf_defaults$mtry
      }
    }
  }

  return(models)
}


#' Validate model specifications
#'
#' @param models List of model specifications
#'
#' @keywords internal
validate_models <- function(models) {

  model_names <- c("Y", "D", "S_e", "S_o")

  for (nm in model_names) {
    if (!nm %in% names(models)) {
      stop("models must contain specifications for Y, D, S_e, and S_o")
    }

    spec <- models[[nm]]

    if (!"model" %in% names(spec)) {
      stop(sprintf("Model specification for %s must include 'model'", nm))
    }

    if (!spec$model %in% c("rf", "logit")) {
      stop(sprintf("Unknown model type for %s: %s", nm, spec$model))
    }

    # Validate RF parameters
    if (spec$model == "rf") {
      if ("num.trees" %in% names(spec) && spec$num.trees < 1) {
        stop(sprintf("num.trees for %s must be positive", nm))
      }
    }
  }

  invisible(TRUE)
}


#' Validate sample split parameters
#'
#' @param train.proportion Training proportion
#'
#' @keywords internal
validate_split_params <- function(train.proportion) {
  if (train.proportion <= 0 || train.proportion >= 1) {
    stop("train.proportion must be between 0 and 1 (exclusive)")
  }
  invisible(TRUE)
}


#' Validate cross-fitting parameters
#'
#' @param nfolds Number of folds
#' @param n Sample size
#'
#' @keywords internal
validate_cv_params <- function(nfolds, n) {
  if (nfolds < 2) {
    stop("nfolds must be at least 2")
  }
  if (nfolds > n) {
    stop("nfolds cannot exceed sample size")
  }
  invisible(TRUE)
}


#' Validate bootstrap parameters
#'
#' @param B Number of bootstrap replications
#' @param clusters Cluster identifiers
#' @param n Sample size
#'
#' @keywords internal
validate_bootstrap_params <- function(B, clusters, n) {
  if (B < 1) {
    stop("B must be positive")
  }
  if (!is.null(clusters) && length(clusters) != n) {
    stop("clusters must have same length as data")
  }
  invisible(TRUE)
}
