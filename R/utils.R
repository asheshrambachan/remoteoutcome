# =============================================================================
# Common Utilities
# =============================================================================

#' Convert binary numeric vector to logical
#'
#' Converts a binary numeric vector (0/1) to a logical vector (TRUE/FALSE).
#' If the input is already logical, it is returned unchanged.
#'
#' @param x A numeric or logical vector. If numeric, must contain only 0 and 1.
#'
#' @return A logical vector.
#'
#' @keywords internal
to_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  } else {
    vals <- unique(x)
    if (all(vals %in% c(0, 1))) {
      return(as.logical(x))
    } else {
      stop("to_logical: expected values 0/1 or TRUE/FALSE; got: ", paste(vals, collapse = ", "))
    }
  }
}


#' Cluster-level bootstrap resampling
#'
#' Performs cluster-level bootstrap resampling by sampling clusters with replacement
#' and returning all observations within the sampled clusters.
#'
#' @param data A data frame containing the data to resample.
#' @param cluster_var Character string naming the cluster variable in \code{data}.
#'   Default is "clusters".
#'
#' @return A data frame with the same structure as \code{data}, containing the
#'   resampled observations.
#'
#' @keywords internal
cluster_sample <- function(data, cluster_var = "clusters") {
  ids <- unique(data[[cluster_var]])
  sampled <- sample(ids, length(ids), replace = TRUE)
  data[data[[cluster_var]] %in% sampled, , drop = FALSE]
}


#' Generate synthetic RSV data
#'
#' Generates synthetic data for testing and simulating RSV estimation methods.
#' The data generating process creates treatment assignment, outcomes, and
#' remotely sensed variables under the RSV framework.
#'
#' @param n Sample size.
#' @param tau Treatment effect parameter.
#' @param X Covariate matrix.
#' @param D Treatment indicator vector.
#' @param Y Outcome vector.
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{Covariate matrix}
#'   \item{D}{Treatment indicator (0/1)}
#'   \item{Y}{Outcome variable (NA for experimental sample)}
#'   \item{Se}{Logical indicator for experimental sample}
#'   \item{So}{Logical indicator for observational sample}
#' }
#'
#' @export
#'
#' @examples
#' # Generate synthetic data
#' n <- 1000
#' X <- matrix(rnorm(n * 5), ncol = 5)
#' D <- rbinom(n, 1, 0.5)
#' Y <- rbinom(n, 1, 0.3)
#' data <- generate_rsv_data(n, tau = 0.1, X, D, Y)
generate_rsv_data <- function(n, tau, X, D, Y) {
  Y1_id <- which(Y == 1)
  Y0_id <- which(Y == 0)
  p_Y1_D1 <- mean(Y[D == 1], na.rm = TRUE) + tau
  p_Y1_D0 <- mean(Y[D == 0], na.rm = TRUE)

  # Draw new D and outcomes
  D <- rbinom(n, 1, prob = 0.5)
  Y_D1 <- rbinom(n, 1, prob = p_Y1_D1)
  Y_D0 <- rbinom(n, 1, prob = p_Y1_D0)
  X_Y1_id <- sample(Y1_id, size = n, replace = TRUE)
  X_Y0_id <- sample(Y0_id, size = n, replace = TRUE)

  # Define Y
  Y <- D * Y_D1 + (1 - D) * Y_D0

  # Define X
  X_Y1 <- X[X_Y1_id, , drop = FALSE]
  X_Y0 <- X[X_Y0_id, , drop = FALSE]
  X <- Y * X_Y1 + (1 - Y) * X_Y0

  # Define Se and So
  Se <- rep(TRUE, n)  # S = e for D = 0 and D = 1
  So <- D != 1        # S = o for D = 0
  Y <- ifelse(So == TRUE, Y, NA)

  list(X = X, D = D, Y = Y, Se = Se, So = So)
}


#' Cross-validation wrapper for RSV estimation
#'
#' Performs k-fold cross-validation for any estimation function that accepts
#' train and test data splits.
#'
#' @param fun Function to apply. Must accept arguments with \code{_train} and
#'   \code{_test} suffixes.
#' @param nfold Number of folds for cross-validation. Default is 5.
#' @param ... Additional arguments passed to \code{fun}. Must include at least
#'   a \code{D} argument to determine sample size.
#'
#' @return A list containing:
#' \describe{
#'   \item{coef_cv}{Cross-validated coefficient estimate}
#'   \item{out}{List of results from each fold}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with RSV estimator
#' cv_results <- cv_rsv(rsv_estimate, nfold = 5,
#'                      X = X, D = D, Y = Y, Se = Se, So = So)
#' }
cv_rsv <- function(fun, nfold = 5, ...) {
  args <- list(...)

  if (!"D" %in% names(args)) stop("cv_rsv: 'D' required")

  n <- length(args$D)
  folds <- sample.int(nfold, n, replace = TRUE)
  out <- vector("list", nfold)

  for (k in seq_len(nfold)) {
    train_id <- folds != k
    test_id <- folds == k

    train_args <- lapply(args, function(el) {
      if (is.matrix(el) || is.data.frame(el)) {
        el[train_id, , drop = FALSE]
      } else {
        el[train_id]
      }
    })
    test_args <- lapply(args, function(el) {
      if (is.matrix(el) || is.data.frame(el)) {
        el[test_id, , drop = FALSE]
      } else {
        el[test_id]
      }
    })

    names(train_args) <- paste0(names(train_args), "_train")
    names(test_args)  <- paste0(names(test_args),  "_test")

    all_args <- c(train_args, test_args)
    call_args <- all_args[names(all_args) %in% names(formals(fun))]
    out[[k]] <- do.call(fun, call_args)
  }

  coef_cv <- mean(sapply(out, function(x) x$coef), na.rm = TRUE)
  list(coef_cv = coef_cv, out = out)
}
