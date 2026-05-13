# =============================================================================
# Bootstrap Standard Errors for RSV Estimates
# =============================================================================

#' Add Bootstrap Standard Errors to RSV Estimate
#'
#' Computes bootstrap standard errors for RSV treatment effect estimates.
#' Works with objects from \code{\link{rsv}}, \code{\link{rsv_split}},
#' \code{\link{rsv_fitted}}, and \code{\link{cv.rsv}}.
#'
#' @param object An object of class \code{"rsv"} or \code{"cv.rsv"}.
#' @param method SE method: \code{"influence"} (analytical, cv.rsv only),
#'   \code{"score_bootstrap"}, or \code{"bootstrap"}.
#' @param B Number of bootstrap replications (default 1000).
#' @param clusters Cluster identifiers for clustered bootstrap (optional).
#'   If \code{NULL}, performs individual-level bootstrap.
#' @param seed Random seed for reproducibility (default \code{NULL}).
#' @param num.threads Number of threads for parallel bootstrap (default 1).
#' @param ... Further arguments (unused; present for S3 method compatibility).
#'
#' @return The input object with added standard error components:
#' \describe{
#'   \item{se}{Standard error of the RSV treatment effect}
#'   \item{relevance_se}{Standard error of the RSV relevance \eqn{E_n[H \cdot \delta_o]}}
#'   \item{se_naive}{Standard error of the naive treatment effect}
#'   \item{relevance_naive_se}{Standard error of the naive relevance}
#'   \item{vcov}{Variance-covariance matrix}
#'   \item{bootstrap}{List with bootstrap information (B, n_failed, seed, tau_draws, seed_draws)}
#' }
#'
#' @details
#' For \code{rsv} objects (from \code{\link{rsv}}, \code{\link{rsv_split}},
#' \code{\link{rsv_fitted}}), the bootstrap resamples observations and
#' recomputes the RSV estimate on each bootstrap sample using the original
#' predictions.
#'
#' For \code{cv.rsv} objects (from \code{\link{cv.rsv}}), the bootstrap
#' resamples observations and recomputes the entire cross-fitting procedure,
#' which is computationally expensive. Consider using fewer replications
#' (e.g., B = 500) for cross-fitted estimates.
#'
#' If \code{clusters} is provided, performs clustered bootstrap by resampling
#' clusters with replacement.
#'
#' @seealso \code{\link{rsv}}, \code{\link{rsv_split}}, \code{\link{rsv_fitted}},
#'   \code{\link{cv.rsv}}
#'
#' @examples
#' \dontrun{
#' # Fit RSV estimate
#' fit <- rsv(Y, D, S_e, S_o, R,
#'            models = list(
#'              Y = list(model = "rf"),
#'              D = list(model = "rf"),
#'              S_e = list(model = "rf"),
#'              S_o = list(model = "rf")
#'            ))
#'
#' # Add bootstrap SE
#' fit <- add_se(fit, B = 1000, num.threads = 8)
#'
#' # With clustered bootstrap
#' fit <- add_se(fit, B = 1000, clusters = village_id, num.threads = 8)
#'
#' # For cross-fitting (expensive!)
#' cv_fit <- cv.rsv(Y, D, S_e, S_o, R,
#'                  models = list(...),
#'                  nfolds = 5)
#' cv_fit <- add_se(cv_fit, B = 500, clusters = village_id, num.threads = 8)
#' }
#'
#' @export
add_se <- function(object, method = c("influence", "score_bootstrap", "bootstrap"),
                   B = 1000, clusters = NULL, seed = NULL, num.threads = 1, ...) {
  UseMethod("add_se")
}


#' @export
add_se.rsv <- function(object, method = c("score_bootstrap"),
                       B = 1000, clusters = NULL, seed = NULL, num.threads = 1, ...) {

  method <- match.arg(method)

  # Store the add_se call
  se_call <- match.call()

  # Extract stored data
  df <- object$data
  theta_init <- object$theta_init
  y_levels <- object$y_levels
  sigma2.lower <- object$sigma2.lower
  sigma2.quantile <- object$sigma2.quantile

  # Validate
  if (is.null(df)) {
    stop("Object does not contain data needed for bootstrap.\n",
         "This may be an old rsv object. Please refit using the updated functions.")
  }

  # Run bootstrap
  boot_result <- .score_bootstrap_rsv(
    df = df,
    theta_init = theta_init,
    y_levels = y_levels,
    sigma2.lower = sigma2.lower,
    sigma2.quantile = sigma2.quantile,
    B = B,
    clusters = clusters,
    seed = seed,
    num.threads = num.threads
  )

  # Add results to object
  object$se                  <- boot_result$se
  object$relevance_se        <- boot_result$relevance_se
  object$se_naive            <- boot_result$se_naive
  object$relevance_naive_se  <- boot_result$relevance_naive_se
  object$vcov                <- boot_result$vcov
  object$bootstrap           <- boot_result$bootstrap
  object$se.method           <- method
  object$se.call             <- se_call

  return(object)
}


#' @export
add_se.cv.rsv <- function(object, method = c("influence", "bootstrap", "score_bootstrap"),
                          B = 1000, clusters = NULL, seed = NULL, num.threads = 1, ...) {

  method <- match.arg(method)

  # Store the add_se call
  se_call <- match.call()

  # Extract stored data
  Y <- object$Y
  D <- object$D
  S_e <- object$S_e
  S_o <- object$S_o
  R <- object$R
  X <- object$X
  X_cols <- object$X_cols
  y_levels <- object$y_levels
  yK <- object$yK
  models <- object$models
  nfolds <- object$nfolds
  fold_ids <- object$fold_ids
  sigma2.lower <- object$sigma2.lower
  sigma2.quantile <- object$sigma2.quantile

  # Validate
  if (is.null(Y) || is.null(D) || is.null(R)) {
    stop("Object does not contain data needed for SE computation.\n",
         "This may be an old cv.rsv object. Please refit using the updated functions.")
  }

  if (method == "influence") {
    # --- Influence function SE (fast, analytical) ---
    pooled_predictions <- object$pooled_predictions
    if (is.null(pooled_predictions)) {
      stop("Object does not contain pooled_predictions needed for influence function SE.\n",
           "Please refit using the latest cv.rsv().")
    }

    result <- .influence_cv_rsv(
      Y = Y, D = D, S_e = S_e, S_o = S_o,
      pooled_predictions = pooled_predictions,
      X = X, X_cols = X_cols,
      y_levels = y_levels, yK = yK,
      sigma2.lower = sigma2.lower,
      sigma2.quantile = sigma2.quantile
    )

    object$se                 <- result$se
    object$vcov               <- result$vcov
    object$influence          <- result$influence
    object$cond_J             <- result$influence$cond_J
    object$se_naive           <- result$se_naive
    object$relevance_se       <- result$relevance_se
    object$relevance_naive_se <- result$relevance_naive_se

  } else if (method == "score_bootstrap") {
    # --- Score bootstrap (fast: resample predictions, no model refitting) ---
    pooled_predictions <- object$pooled_predictions
    if (is.null(pooled_predictions)) {
      stop("Object does not contain pooled_predictions needed for score bootstrap.\n",
           "Please refit using the latest cv.rsv().")
    }

    boot_result <- .score_bootstrap_cv_rsv(
      Y = Y, D = D, S_e = S_e, S_o = S_o,
      pooled_predictions = pooled_predictions,
      X = X, X_cols = X_cols,
      y_levels = y_levels, yK = yK,
      sigma2.lower = sigma2.lower,
      sigma2.quantile = sigma2.quantile,
      B = B, clusters = clusters, seed = seed,
      num.threads = num.threads
    )

    object$se                  <- boot_result$se
    object$relevance_se        <- boot_result$relevance_se
    object$se_naive            <- boot_result$se_naive
    object$relevance_naive_se  <- boot_result$relevance_naive_se
    object$vcov                <- boot_result$vcov
    object$bootstrap           <- boot_result$bootstrap
    object$cond_J_median       <- boot_result$bootstrap$kappa_J_median
    object$cond_J_max          <- boot_result$bootstrap$kappa_J_max

  } else {
    # --- Bootstrap SE (expensive: refits models in every fold × replicate) ---
    boot_result <- .bootstrap_cv_rsv(
      Y = Y, D = D, S_e = S_e, S_o = S_o, R = R,
      X = X, X_cols = X_cols,
      y_levels = y_levels, yK = yK,
      models = models, nfolds = nfolds, fold_ids = fold_ids,
      sigma2.lower = sigma2.lower,
      sigma2.quantile = sigma2.quantile,
      B = B, clusters = clusters, seed = seed,
      num.threads = num.threads
    )

    object$se                  <- boot_result$se
    object$relevance_se        <- boot_result$relevance_se
    object$se_naive            <- boot_result$se_naive
    object$relevance_naive_se  <- boot_result$relevance_naive_se
    object$vcov                <- boot_result$vcov
    object$bootstrap           <- boot_result$bootstrap
  }

  object$se.method <- method
  object$se.call <- se_call
  return(object)
}


# =============================================================================
# Internal Bootstrap Functions
# =============================================================================

#' Bootstrap RSV estimate from combined data frame
#'
#' @param df Data frame with observations (Y, D, S_e, S_o, X) and predictions (pred_*)
#' @param theta_init Initial treatment effect estimate
#' @param sigma2.lower Lower bound for sigma-squared
#' @param sigma2.quantile Quantile for sigma-squared lower bound
#' @param B Number of bootstrap replications
#' @param clusters Cluster identifiers (optional)
#' @param seed Random seed (optional)
#' @param num.threads Number of threads for parallel bootstrap
#'
#' @return List with se, relevance_se, vcov, and bootstrap info
#'
#' @keywords internal
.score_bootstrap_rsv <- function(df, theta_init, y_levels,
                           sigma2.lower, sigma2.quantile,
                           B, clusters, seed, num.threads) {

  n <- nrow(df)
  validate_bootstrap_params(B, clusters, n)

  # Setup clusters
  if (is.null(clusters)) {
    clusters <- seq_len(n)
    unique_clusters <- clusters
  } else {
    unique_clusters <- unique(clusters)
  }

  # Bootstrap function for one replication
  run_one <- function(b) {

    # Set seed if provided
    if (!is.null(seed)) {
      set.seed(seed + b)
    }

    # Cluster bootstrap: resample clusters with replacement
    clusters_b <- sample(unique_clusters, length(unique_clusters), replace = TRUE)

    # Create index mapping for resampled data
    idx_b <- unlist(lapply(clusters_b, function(cl) which(clusters == cl)))

    # Resample combined data
    df_b <- df[idx_b, , drop = FALSE]

    # Compute estimate on bootstrap sample
    tryCatch({
      est_b <- compute_estimate(
        df = df_b,
        theta_init = theta_init,
        y_levels = y_levels,
        sigma2.lower = sigma2.lower,
        sigma2.quantile = sigma2.quantile
      )
      # Naive estimate
      est_naive_b <- compute_estimate(
        df              = df_b,
        theta_init      = theta_init,
        y_levels        = y_levels,
        sigma2.lower    = sigma2.lower,
        sigma2.quantile = sigma2.quantile,
        pred_y_only     = TRUE
      )
      c(coef = unname(est_b$coefficients), denominator = est_b$relevance_ate,
        coef_naive = unname(est_naive_b$coefficients), denominator_naive = est_naive_b$relevance_ate)
    }, error = function(e) {
      c(coef = NA_real_, denominator = NA_real_, coef_naive = NA_real_, denominator_naive = NA_real_)
    })
  }

  # Run bootstrap
  boot_results <- parallel::mclapply(seq_len(B), run_one, mc.cores = num.threads)

  # Extract results
  boot_coefs        <- vapply(boot_results, function(x) x["coef"],             FUN.VALUE = numeric(1))
  boot_denoms       <- vapply(boot_results, function(x) x["denominator"],      FUN.VALUE = numeric(1))
  boot_coefs_naive  <- vapply(boot_results, function(x) x["coef_naive"],       FUN.VALUE = numeric(1))
  boot_denoms_naive <- vapply(boot_results, function(x) x["denominator_naive"],FUN.VALUE = numeric(1))

  # Compute standard errors
  se                 <- stats::sd(boot_coefs,        na.rm = TRUE)
  relevance_se       <- stats::sd(boot_denoms,       na.rm = TRUE)
  se_naive           <- stats::sd(boot_coefs_naive,  na.rm = TRUE)
  relevance_naive_se <- stats::sd(boot_denoms_naive, na.rm = TRUE)

  # Check for failures
  n_failed <- sum(is.na(boot_coefs))
  if (n_failed > 0) {
    warning(sprintf("%d of %d bootstrap replications failed", n_failed, B))
  }

  list(
    se                 = structure(se,       names = "D"),
    relevance_se       = relevance_se,
    se_naive           = structure(se_naive, names = "D"),
    relevance_naive_se = relevance_naive_se,
    vcov = structure(
      matrix(se^2, 1, 1),
      dimnames = list("D", "D")
    ),
    bootstrap = list(
      B        = B,
      n_failed = n_failed,
      clusters = if (is.null(clusters)) NULL else "provided",
      seed      = seed,
      tau_draws = boot_coefs
    )
  )
}


#' Score bootstrap for cross-fitted RSV estimate
#'
#' Fast bootstrap that resamples (obs, prediction) pairs without refitting any
#' models.  Since predictions are fixed from the original cv.rsv fit, fold
#' assignments carry no information and all n observations are resampled
#' together.  No between-fold variance correction is applied.
#'
#' @keywords internal
.score_bootstrap_cv_rsv <- function(Y, D, S_e, S_o,
                                    pooled_predictions,
                                    X, X_cols, y_levels, yK,
                                    sigma2.lower, sigma2.quantile,
                                    B, clusters, seed, num.threads) {

  n <- length(Y)
  validate_bootstrap_params(B, clusters, n)

  # Build pooled df (obs + fixed predictions)
  df_pooled <- cbind(
    data.frame(Y = Y, D = D, S_e = S_e, S_o = S_o),
    pooled_predictions
  )
  df_pooled <- .ensure_numeric_indicators(df_pooled)

  X_cols_use <- if (!is.null(X_cols) && length(X_cols) > 0) X_cols else character(0)
  theta_init <- compute_theta_init(df_pooled, y_levels = y_levels,
                                   yK = yK, X_cols = X_cols_use)

  # Setup clusters
  if (is.null(clusters)) {
    cluster_ids     <- seq_len(n)
    unique_clusters <- cluster_ids
  } else {
    cluster_ids     <- clusters
    unique_clusters <- unique(clusters)
  }

  C      <- length(unique_clusters)
  cl_idx <- match(cluster_ids, unique_clusters)
  cl_sizes <- tabulate(cl_idx, nbins = C)

  run_one <- function(b) {
    seed_b <- if (!is.null(seed)) seed + b else NA_integer_
    if (!is.na(seed_b)) set.seed(seed_b)

    clusters_b <- sample(unique_clusters, C, replace = TRUE)
    idx_b      <- unlist(lapply(clusters_b, function(cl) which(cluster_ids == cl)))
    df_use     <- df_pooled[idx_b, , drop = FALSE]
    w_use      <- NULL

    tryCatch({
      # RSV estimate
      est_b <- compute_estimate(
        df              = df_use,
        theta_init      = theta_init,
        y_levels        = y_levels,
        sigma2.lower    = sigma2.lower,
        sigma2.quantile = sigma2.quantile,
        weights         = w_use
      )
      denom     <- est_b$denominator$denominator
      rel_val   <- est_b$relevance_ate

      # kappa(J) for this resample: one J matrix per X stratum, take the max
      kappa_J_b <- max(vapply(denom, function(J)
        if (is.matrix(J) && nrow(J) > 1 && all(is.finite(J))) kappa(J, exact = TRUE) else NA_real_,
        FUN.VALUE = numeric(1)), na.rm = TRUE)

      # Naive estimate
      est_naive_b <- compute_estimate(
        df              = df_use,
        theta_init      = theta_init,
        y_levels        = y_levels,
        sigma2.lower    = sigma2.lower,
        sigma2.quantile = sigma2.quantile,
        pred_y_only     = TRUE,
        weights         = w_use
      )
      rel_naive_val <- est_naive_b$relevance_ate

      c(seed_b       = seed_b,
        coef         = unname(est_b$coefficients),
        relevance    = rel_val,
        coef_naive   = unname(est_naive_b$coefficients),
        relevance_naive = rel_naive_val,
        kappa_J      = kappa_J_b)
    }, error = function(e) {
      c(seed_b = seed_b, coef = NA_real_, relevance = NA_real_,
        coef_naive = NA_real_, relevance_naive = NA_real_,
        kappa_J = NA_real_)
    })
  }

  boot_results      <- parallel::mclapply(seq_len(B), run_one, mc.cores = num.threads)
  boot_coefs        <- vapply(boot_results, function(x) x["coef"],            FUN.VALUE = numeric(1))
  boot_rels         <- vapply(boot_results, function(x) x["relevance"],       FUN.VALUE = numeric(1))
  boot_coefs_naive  <- vapply(boot_results, function(x) x["coef_naive"],      FUN.VALUE = numeric(1))
  boot_rels_naive   <- vapply(boot_results, function(x) x["relevance_naive"], FUN.VALUE = numeric(1))
  boot_seeds        <- vapply(boot_results, function(x) x["seed_b"],          FUN.VALUE = numeric(1))
  boot_kappas       <- vapply(boot_results, function(x) x["kappa_J"],         FUN.VALUE = numeric(1))

  se                <- stats::sd(boot_coefs,       na.rm = TRUE)
  relevance_se      <- stats::sd(boot_rels,        na.rm = TRUE)
  se_naive          <- stats::sd(boot_coefs_naive, na.rm = TRUE)
  relevance_naive_se <- stats::sd(boot_rels_naive, na.rm = TRUE)

  n_failed <- sum(is.na(boot_coefs))
  if (n_failed > 0)
    warning(sprintf("%d of %d score bootstrap replications failed", n_failed, B))

  list(
    se                 = structure(se,       names = "D"),
    relevance_se       = relevance_se,
    se_naive           = structure(se_naive, names = "D"),
    relevance_naive_se = relevance_naive_se,
    vcov               = structure(matrix(se^2, 1, 1), dimnames = list("D", "D")),
    bootstrap = list(
      B              = B,
      n_failed       = n_failed,
      clusters       = if (is.null(clusters)) NULL else "provided",
      seed           = seed,
      tau_draws      = boot_coefs,
      seed_draws     = boot_seeds,
      kappa_J_draws  = boot_kappas,
      kappa_J_median = median(boot_kappas, na.rm = TRUE),
      kappa_J_max    = max(boot_kappas,    na.rm = TRUE)
    )
  )
}


#' Bootstrap cross-fitted RSV estimate
#'
#' Stratified bootstrap: resamples independently within each fold so that fold
#' sizes are preserved and train/test leakage is prevented.  Models are refitted
#' on each bootstrap replicate.  After collecting test-fold predictions across
#' all folds, a single \code{compute_estimate} call is made on the full pooled
#' bootstrap sample — mirroring the structure of \code{.score_bootstrap_cv_rsv}
#' with the only difference being refitted vs. fixed predictions.
#'
#' @keywords internal
.bootstrap_cv_rsv <- function(Y, D, S_e, S_o, R,
                              X, X_cols, y_levels, yK,
                              models, nfolds, fold_ids,
                              sigma2.lower, sigma2.quantile,
                              B, clusters, seed, num.threads) {

  n <- length(Y)
  validate_bootstrap_params(B, clusters, n)

  X_cols_use <- if (!is.null(X_cols) && length(X_cols) > 0) X_cols else character(0)

  # Setup clusters
  if (is.null(clusters)) {
    cluster_ids     <- seq_len(n)
    unique_clusters <- cluster_ids
  } else {
    cluster_ids     <- clusters
    unique_clusters <- unique(clusters)
  }

  # Bootstrap function for one replication
  run_one <- function(b) {
    if (!is.null(seed)) set.seed(seed + b)

    # Stratified within-fold resampling: resample independently inside each
    # fold so that every bootstrap sample has exactly the same fold sizes as
    # the original, preventing train/test leakage.
    idx_b <- unlist(lapply(seq_len(nfolds), function(k) {
      fold_k_obs <- which(fold_ids == k)
      if (is.null(clusters)) {
        sample(fold_k_obs, length(fold_k_obs), replace = TRUE)
      } else {
        cls_in_fold_k <- unique(cluster_ids[fold_k_obs])
        cls_b_k <- sample(cls_in_fold_k, length(cls_in_fold_k), replace = TRUE)
        unlist(lapply(cls_b_k, function(cl) which(cluster_ids == cl)))
      }
    }))

    Y_b        <- Y[idx_b]
    D_b        <- D[idx_b]
    S_e_b      <- S_e[idx_b]
    S_o_b      <- S_o[idx_b]
    R_b        <- R[idx_b, , drop = FALSE]
    X_b        <- if (!is.null(X)) X[idx_b, , drop = FALSE] else NULL
    fold_ids_b <- fold_ids[idx_b]

    tryCatch({

      # Refit models per fold; collect test-fold (obs + predictions) rows
      fold_dfs <- vector("list", nfolds)

      for (k in seq_len(nfolds)) {
        train_idx <- (fold_ids_b != k)
        test_idx  <- (fold_ids_b == k)

        fit_result <- fit_all_models(
          R_train = R_b[train_idx, , drop = FALSE],
          X_train = if (!is.null(X_b)) X_b[train_idx, , drop = FALSE] else NULL,
          Y = Y_b[train_idx],
          D = D_b[train_idx],
          S_e = S_e_b[train_idx],
          S_o = S_o_b[train_idx],
          R_test = R_b[test_idx, , drop = FALSE],
          X_test = if (!is.null(X_b)) X_b[test_idx, , drop = FALSE] else NULL,
          X_cols = X_cols,
          y_levels = y_levels,
          yK = yK,
          models = models,
          seed = NULL,
          num.threads = 1
        )

        fold_dfs[[k]] <- cbind(
          data.frame(
            Y   = Y_b[test_idx],
            D   = D_b[test_idx],
            S_e = S_e_b[test_idx],
            S_o = S_o_b[test_idx]
          ),
          fit_result$predictions
        )
      }

      # Pool all test-fold predictions and compute one estimate on all n obs,
      # mirroring .score_bootstrap_cv_rsv's single compute_estimate call
      df_b         <- do.call(rbind, fold_dfs)
      theta_init_b <- compute_theta_init(df_b, y_levels = y_levels,
                                         yK = yK, X_cols = X_cols_use)

      # RSV estimate
      est_b <- compute_estimate(
        df              = df_b,
        theta_init      = theta_init_b,
        y_levels        = y_levels,
        sigma2.lower    = sigma2.lower,
        sigma2.quantile = sigma2.quantile
      )
      rel_val <- est_b$relevance_ate

      # Naive estimate
      est_naive_b <- compute_estimate(
        df              = df_b,
        theta_init      = theta_init_b,
        y_levels        = y_levels,
        sigma2.lower    = sigma2.lower,
        sigma2.quantile = sigma2.quantile,
        pred_y_only     = TRUE
      )
      rel_naive_val <- est_naive_b$relevance_ate

      c(seed_b         = if (!is.null(seed)) seed + b else NA_real_,
        coef           = unname(est_b$coefficients),
        relevance      = rel_val,
        coef_naive     = unname(est_naive_b$coefficients),
        relevance_naive = rel_naive_val)

    }, error = function(e) {
      c(seed_b = if (!is.null(seed)) seed + b else NA_real_,
        coef = NA_real_, relevance = NA_real_,
        coef_naive = NA_real_, relevance_naive = NA_real_)
    })
  }

  boot_results      <- parallel::mclapply(seq_len(B), run_one, mc.cores = num.threads)
  boot_coefs        <- vapply(boot_results, function(x) x["coef"],            FUN.VALUE = numeric(1))
  boot_rels         <- vapply(boot_results, function(x) x["relevance"],       FUN.VALUE = numeric(1))
  boot_coefs_naive  <- vapply(boot_results, function(x) x["coef_naive"],      FUN.VALUE = numeric(1))
  boot_rels_naive   <- vapply(boot_results, function(x) x["relevance_naive"], FUN.VALUE = numeric(1))
  boot_seeds        <- vapply(boot_results, function(x) x["seed_b"],          FUN.VALUE = numeric(1))

  se                 <- stats::sd(boot_coefs,       na.rm = TRUE)
  relevance_se       <- stats::sd(boot_rels,        na.rm = TRUE)
  se_naive           <- stats::sd(boot_coefs_naive, na.rm = TRUE)
  relevance_naive_se <- stats::sd(boot_rels_naive,  na.rm = TRUE)

  n_failed <- sum(is.na(boot_coefs))
  if (n_failed > 0)
    warning(sprintf("%d of %d bootstrap replications failed", n_failed, B))

  list(
    se                 = structure(se,       names = "D"),
    relevance_se       = relevance_se,
    se_naive           = structure(se_naive, names = "D"),
    relevance_naive_se = relevance_naive_se,
    vcov               = structure(matrix(se^2, 1, 1), dimnames = list("D", "D")),
    bootstrap = list(
      B          = B,
      n_failed   = n_failed,
      clusters   = if (is.null(clusters)) NULL else "provided",
      seed       = seed,
      tau_draws  = boot_coefs,
      seed_draws = boot_seeds
    )
  )
}


#' Influence function SE for cross-fitted RSV estimator
#'
#' Computes analytical standard errors via the asymptotic linear representation
#' of the cross-fitted RSV estimator (Proposition C.4). Implements the full
#' delta-method approach for general discrete outcomes (K >= 2), following
#' misc/v2/crossfit_standard_errors.tex.
#'
#' Structure (matching tex notation):
#'   For each covariate cell x (or the full sample when X_cols is empty):
#'     1. Assemble \eqn{\phi_i(x)}: raw score before applying \eqn{A(x) = J(x)^{-1}}
#'          \eqn{\phi_i = H_i r_i} + delta-method corrections for p_e1, p_e0, p_oj, p_oK
#'     2. Apply \eqn{A(x) = J(x)^{-1}} once: \eqn{\psi_{\theta,i}(x) = J(x)^{-1} \phi_i(x)}
#'     3. Project to scalar ATE scale: \eqn{\psi_{\tau,i}(x) = \lambda' \psi_{\theta,i}(x)}
#'
#'   No-covariate case (X_cols empty):
#'     \eqn{\hat{\psi}_i = \lambda' J^{-1} \phi_i}
#'
#'   Covariate case (X_cols non-empty), per tex eq:
#'     \eqn{\hat{\psi}_i = (\pi(x_i)/P(X=x_i)) \lambda' A(x_i) \phi_i(x_i)}
#'               \eqn{+ P(S=e)^{-1} \mathbf{1}\{S_i=e\} (\tau(x_i) - \bar{\tau})}
#'     where \eqn{\bar{\tau} = \sum_x \pi(x) \tau(x)} (experimental average ATE)
#'
#' SE = sqrt( mean(psi_i^2) / n )
#'
#' @param Y Outcome vector
#' @param D Treatment vector
#' @param S_e Experimental sample indicator
#' @param S_o Observational sample indicator
#' @param pooled_predictions Data frame of cross-fitted predictions (n rows, original order)
#' @param X Discrete covariates (optional)
#' @param X_cols Column names for X
#' @param y_levels Levels of Y (length K)
#' @param yK Reference category
#' @param sigma2.lower Lower bound for sigma2
#' @param sigma2.quantile Quantile for sigma2 lower bound
#'
#' @return List with coef, se, se_naive, relevance_se, relevance_naive_se, vcov, influence
#'
#' @keywords internal
.influence_cv_rsv <- function(Y, D, S_e, S_o,
                              pooled_predictions,
                              X, X_cols, y_levels, yK,
                              sigma2.lower, sigma2.quantile) {

  n <- length(Y)
  if (is.null(yK)) yK <- y_levels[1]

  yj_levels <- setdiff(y_levels, yK)
  M         <- length(yj_levels)
  lambda    <- as.numeric(yj_levels) - as.numeric(yK)   # M-vector: yj - yK

  # -----------------------------------------------------------------------
  # 1) Build pooled df, compute RSV weights H and naive weights H_naive
  #    H[i,j]       = delta_o_pred[i,j] / sigma2[i]
  #    H_naive[i,j] = pred_Y_j[i]
  # -----------------------------------------------------------------------
  df_pooled <- cbind(
    data.frame(Y = Y, D = D, S_e = S_e, S_o = S_o),
    pooled_predictions
  )
  X_cols_use <- if (!is.null(X_cols) && length(X_cols) > 0) X_cols else character(0)
  df_pooled <- .ensure_numeric_indicators(df_pooled)

  theta_init <- compute_theta_init(df_pooled, y_levels = y_levels,
                                   yK = yK, X_cols = X_cols_use)
  H_mat <- as.matrix(
    compute_H(df_pooled, theta_init, y_levels = y_levels,
              sigma2.lower = sigma2.lower,
              sigma2.quantile = sigma2.quantile)
  )  # n x M

  H_naive_mat <- as.matrix(
    compute_H(df_pooled, theta_init, y_levels = y_levels,
              sigma2.lower = sigma2.lower,
              sigma2.quantile = sigma2.quantile,
              pred_y_only = TRUE)
  )  # n x M

  # -----------------------------------------------------------------------
  # 2) Observed indicator variables  (NA & FALSE = FALSE, so NA-safe)
  # -----------------------------------------------------------------------
  ind_D1Se  <- as.numeric(D == 1 & S_e == 1)   # 1{D=1, S=e}
  ind_D0Se  <- as.numeric(D == 0 & S_e == 1)   # 1{D=0, S=e}
  ind_Yj_So <- matrix(0, n, M)                 # col j = 1{Y=yj, S=o}
  for (j in seq_len(M)) ind_Yj_So[, j] <- as.numeric(Y == yj_levels[j] & S_o == 1)
  ind_YK_So <- as.numeric(Y == yK & S_o == 1)  # 1{Y=yK, S=o}

  # -----------------------------------------------------------------------
  # 3) Per-cell influence for RSV and naive (tau and relevance)
  #
  #    phi_i(x) = H_i r_i + delta-method corrections for p_e1, p_e0, p_oj, p_oK
  #    psi_theta_i(x) = J(x)^{-1} phi_i(x)
  #
  #    For relevance rho = 1/(lambda' J^{-1} lambda), g = J^{-1} lambda:
  #    psi_rho_i(x) = rho^2 (H_i' g)(delta_o_i' g) - rho
  #                 + sum_j [-rho^2 g_j (v_oj_j' g)] (Yj_So_i - p_oj)
  #                 + rho^2 sum(g) (g' v_oK) (YK_So_i - p_oK)
  #
  #    Naive quantities: same formulas with H_naive replacing H.
  # -----------------------------------------------------------------------
  cell_influence <- function(idx) {
    n_x      <- length(idx)
    H_x      <- H_mat[idx, , drop = FALSE]        # n_x x M (RSV)
    H_nv_x   <- H_naive_mat[idx, , drop = FALSE]  # n_x x M (naive)
    D1Se_x   <- ind_D1Se[idx]
    D0Se_x   <- ind_D0Se[idx]
    Yj_x     <- ind_Yj_So[idx, , drop = FALSE]    # n_x x M
    YK_x     <- ind_YK_So[idx]

    p_e1 <- mean(D1Se_x)
    p_e0 <- mean(D0Se_x)
    p_oj <- colMeans(Yj_x)
    p_oK <- mean(YK_x)

    delta_e <- D1Se_x / p_e1 - D0Se_x / p_e0
    delta_o <- sweep(Yj_x, 2, p_oj, "/") - YK_x * (1 / p_oK)    # n_x x M

    # ---- RSV: J, theta, phi, psi_theta, psi_rho ----
    J_x     <- crossprod(H_x, delta_o) / n_x
    a_x     <- drop(crossprod(H_x, delta_e)) / n_x
    kappa_J <- if (all(is.finite(J_x))) kappa(J_x, exact = TRUE) else NA_real_
    J_x_inv <- solve(J_x)
    theta_x <- drop(J_x_inv %*% a_x)
    r_x     <- delta_e - drop(delta_o %*% theta_x)

    g_rsv   <- drop(J_x_inv %*% lambda)
    rho_lq  <- drop(t(lambda) %*% J_x_inv %*% lambda)
    rho_rsv <- if (is.finite(rho_lq) && rho_lq != 0) 1 / rho_lq else NA_real_

    v_e1 <- colMeans(H_x * D1Se_x) / p_e1^2
    v_e0 <- colMeans(H_x * D0Se_x) / p_e0^2
    v_oK <- colMeans(H_x * YK_x)   / p_oK^2

    phi     <- H_x * r_x
    phi     <- phi - (D1Se_x - p_e1) %o% v_e1
    phi     <- phi + (D0Se_x - p_e0) %o% v_e0
    psi_rho <- if (!is.na(rho_rsv)) {
      rho_rsv^2 * drop(H_x %*% g_rsv) * drop(delta_o %*% g_rsv) - rho_rsv
    } else rep(NA_real_, n_x)

    for (j in seq_len(M)) {
      v_oj_j <- colMeans(H_x * Yj_x[, j]) / p_oj[j]^2
      phi    <- phi + theta_x[j] * (Yj_x[, j] - p_oj[j]) %o% v_oj_j
      if (!is.na(rho_rsv))
        psi_rho <- psi_rho -
          rho_rsv^2 * drop(g_rsv %*% v_oj_j) * g_rsv[j] * (Yj_x[, j] - p_oj[j])
    }
    phi <- phi - sum(theta_x) * (YK_x - p_oK) %o% v_oK
    if (!is.na(rho_rsv))
      psi_rho <- psi_rho +
        rho_rsv^2 * drop(g_rsv %*% v_oK) * sum(g_rsv) * (YK_x - p_oK)

    psi_theta <- t(J_x_inv %*% t(phi))   # n_x x M

    # ---- Naive: J_nv, theta_nv, phi_nv, psi_theta_nv, psi_rho_nv ----
    J_nv_x     <- crossprod(H_nv_x, delta_o) / n_x
    a_nv_x     <- drop(crossprod(H_nv_x, delta_e)) / n_x
    J_nv_inv   <- solve(J_nv_x)
    theta_nv_x <- drop(J_nv_inv %*% a_nv_x)
    r_nv_x     <- delta_e - drop(delta_o %*% theta_nv_x)

    g_nv      <- drop(J_nv_inv %*% lambda)
    rho_nv_lq <- drop(t(lambda) %*% J_nv_inv %*% lambda)
    rho_nv    <- if (is.finite(rho_nv_lq) && rho_nv_lq != 0) 1 / rho_nv_lq else NA_real_

    v_e1_nv <- colMeans(H_nv_x * D1Se_x) / p_e1^2
    v_e0_nv <- colMeans(H_nv_x * D0Se_x) / p_e0^2
    v_oK_nv <- colMeans(H_nv_x * YK_x)   / p_oK^2

    phi_nv     <- H_nv_x * r_nv_x
    phi_nv     <- phi_nv - (D1Se_x - p_e1) %o% v_e1_nv
    phi_nv     <- phi_nv + (D0Se_x - p_e0) %o% v_e0_nv
    psi_rho_nv <- if (!is.na(rho_nv)) {
      rho_nv^2 * drop(H_nv_x %*% g_nv) * drop(delta_o %*% g_nv) - rho_nv
    } else rep(NA_real_, n_x)

    for (j in seq_len(M)) {
      v_oj_j_nv <- colMeans(H_nv_x * Yj_x[, j]) / p_oj[j]^2
      phi_nv    <- phi_nv + theta_nv_x[j] * (Yj_x[, j] - p_oj[j]) %o% v_oj_j_nv
      if (!is.na(rho_nv))
        psi_rho_nv <- psi_rho_nv -
          rho_nv^2 * drop(g_nv %*% v_oj_j_nv) * g_nv[j] * (Yj_x[, j] - p_oj[j])
    }
    phi_nv <- phi_nv - sum(theta_nv_x) * (YK_x - p_oK) %o% v_oK_nv
    if (!is.na(rho_nv))
      psi_rho_nv <- psi_rho_nv +
        rho_nv^2 * drop(g_nv %*% v_oK_nv) * sum(g_nv) * (YK_x - p_oK)

    psi_theta_nv <- t(J_nv_inv %*% t(phi_nv))   # n_x x M

    list(psi_theta   = psi_theta,
         psi_rho     = psi_rho,
         psi_tau_nv  = drop(psi_theta_nv %*% lambda),
         psi_rho_nv  = psi_rho_nv,
         theta_x     = theta_x,
         tau_x       = sum(lambda * theta_x),
         rho_x       = rho_rsv,
         tau_nv_x    = sum(lambda * theta_nv_x),
         rho_nv_x    = rho_nv,
         J_x         = J_x,
         J_x_inv     = J_x_inv,
         kappa_J     = kappa_J)
  }

  # -----------------------------------------------------------------------
  # 4a) No-covariate case: single cell over all observations
  # -----------------------------------------------------------------------
  if (length(X_cols_use) == 0) {

    res           <- cell_influence(seq_len(n))
    psi_theta_mat <- res$psi_theta
    theta_hat     <- res$theta_x
    theta0_pooled <- res$tau_x

    psi_tau0   <- drop(psi_theta_mat %*% lambda)
    psi_rho0   <- res$psi_rho
    psi_tau_nv <- res$psi_tau_nv
    psi_rho_nv <- res$psi_rho_nv

    inf_out <- list(
      psi_tau0  = psi_tau0,
      psi_theta = psi_theta_mat,
      theta_hat = theta_hat,
      J         = res$J_x,
      J_inv     = res$J_x_inv,
      lambda    = lambda,
      cond_J    = res$kappa_J
    )

  # -----------------------------------------------------------------------
  # 4b) Covariate case: loop over cells x in X
  #
  #   psi_tau_i = [pi(x_i)/P(X=x_i)] lambda' J(x_i)^{-1} phi_i(x_i)
  #             + P(S=e)^{-1} 1{S_i=e} (tau(x_i) - tau_bar)
  #
  #   psi_rho_i = [pi(x_i)/P(X=x_i)] psi_rho_i(x_i)
  #             + P(S=e)^{-1} 1{S_i=e} (rho(x_i) - rho_bar)
  #
  #   (Same structure for naive quantities.)
  # -----------------------------------------------------------------------
  } else {

    X_df         <- as.data.frame(X)[, X_cols_use, drop = FALSE]
    cell_labels  <- do.call(paste, c(X_df, sep = "\u001f"))   # unit-sep avoids collisions
    unique_cells <- unique(cell_labels)

    P_Se          <- mean(S_e == 1)
    psi_tau0      <- numeric(n)
    psi_rho0      <- numeric(n)
    psi_tau_nv    <- numeric(n)
    psi_rho_nv    <- numeric(n)
    tau_by_obs    <- numeric(n)
    rho_by_obs    <- numeric(n)
    tau_nv_by_obs <- numeric(n)
    rho_nv_by_obs <- numeric(n)
    cell_kappas   <- numeric(length(unique_cells))

    for (ci in seq_along(unique_cells)) {
      cx    <- unique_cells[ci]
      idx_x <- which(cell_labels == cx)
      n_x   <- length(idx_x)
      res   <- cell_influence(idx_x)

      tau_by_obs[idx_x]    <- res$tau_x
      rho_by_obs[idx_x]    <- res$rho_x
      tau_nv_by_obs[idx_x] <- res$tau_nv_x
      rho_nv_by_obs[idx_x] <- res$rho_nv_x
      cell_kappas[ci]      <- res$kappa_J

      n_xe  <- sum(S_e[idx_x] == 1)
      pi_x  <- (n_xe / n) / P_Se
      P_X_x <- n_x / n
      w_x   <- pi_x / P_X_x

      psi_tau0[idx_x]   <- psi_tau0[idx_x]   + w_x * drop(res$psi_theta %*% lambda)
      psi_rho0[idx_x]   <- psi_rho0[idx_x]   + w_x * res$psi_rho
      psi_tau_nv[idx_x] <- psi_tau_nv[idx_x] + w_x * res$psi_tau_nv
      psi_rho_nv[idx_x] <- psi_rho_nv[idx_x] + w_x * res$psi_rho_nv
    }

    Se_ind        <- as.numeric(S_e == 1)
    tau_bar       <- sum(tau_by_obs[S_e == 1])    / sum(S_e == 1)
    rho_bar       <- sum(rho_by_obs[S_e == 1])    / sum(S_e == 1)
    tau_nv_bar    <- sum(tau_nv_by_obs[S_e == 1]) / sum(S_e == 1)
    rho_nv_bar    <- sum(rho_nv_by_obs[S_e == 1]) / sum(S_e == 1)
    theta0_pooled <- tau_bar

    psi_tau0   <- psi_tau0   + Se_ind * (tau_by_obs    - tau_bar)    / P_Se
    psi_rho0   <- psi_rho0   + Se_ind * (rho_by_obs    - rho_bar)    / P_Se
    psi_tau_nv <- psi_tau_nv + Se_ind * (tau_nv_by_obs - tau_nv_bar) / P_Se
    psi_rho_nv <- psi_rho_nv + Se_ind * (rho_nv_by_obs - rho_nv_bar) / P_Se

    max_kappa_J <- if (all(is.na(cell_kappas))) NA_real_ else max(cell_kappas, na.rm = TRUE)
    inf_out <- list(
      psi_tau0   = psi_tau0,
      tau_by_obs = tau_by_obs,
      tau_bar    = tau_bar,
      lambda     = lambda,
      cond_J     = max_kappa_J
    )
  }

  # -----------------------------------------------------------------------
  # 5) SE = sqrt( mean(psi_i^2) / n ) for each quantity
  # -----------------------------------------------------------------------
  se_val                 <- sqrt(mean(psi_tau0^2)                / n)
  se_naive_val           <- sqrt(mean(psi_tau_nv^2)              / n)
  relevance_se_val       <- sqrt(mean(psi_rho0^2,  na.rm = TRUE) / n)
  relevance_naive_se_val <- sqrt(mean(psi_rho_nv^2, na.rm = TRUE) / n)

  list(
    coef               = structure(theta0_pooled, names = "D"),
    se                 = structure(se_val,       names = "D"),
    se_naive           = structure(se_naive_val, names = "D"),
    relevance_se       = relevance_se_val,
    relevance_naive_se = relevance_naive_se_val,
    vcov = structure(
      matrix(se_val^2, 1, 1),
      dimnames = list("D", "D")
    ),
    influence = inf_out
  )
}
