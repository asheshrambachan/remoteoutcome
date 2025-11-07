# =============================================================================
# RSV Helper Functions (Binary Outcome)
# Following main text Algorithm 1 with Remark 3 generalization
# =============================================================================

#' Fit predictions using random forest
#'
#' Fits random forest models to predict Y(R), D(R), S_e(R), and S_o(R) using
#' the remotely sensed variable R as the only predictor.
#'
#' @param R Remotely sensed variable.
#' @param Y Outcome variable (binary, NA for experimental-only sample).
#' @param D Treatment indicator (binary, NA for observational-only sample).
#' @param S_e Experimental sample indicator (0 or 1).
#' @param S_o Observational sample indicator (0 or 1).
#' @param R_pred Vector of R values on which to make predictions (defaults to R).
#' @param ml_params List of parameters for random forest. Should include:
#'   \describe{
#'     \item{ntree}{Number of trees}
#'     \item{classwt_Y}{Class weights for \code{pred_Y} model; default \code{c(10, 1)})}
#'     \item{seed}{User specified seed passed to each \code{ranger} function for reproducibility; default \code{NULL}}
#'     \item{cores}{Number of cores used by `ranger` for parallel training; default \code{1}}
#'   }
#'
#' @return A list containing: :
#'   \describe{
#'     \item{theta_init}{Initial estimate of the treatment effect}
#'     \item{predictions}{A data frame with columns: \code{pred_Y}}, \code{pred_D}, \code{pred_S_e}, \code{pred_S_o}
#'   }
#'
#' @keywords internal
fit_predictions_rf <- function(R, Y, D, S_e, S_o, R_pred = NULL, ml_params = list()) {
  if (is.null(R_pred)) R_pred <- R
  
  if (class(S_e) == "logical")
    S_e <- as.numeric(S_e)
  
  if (class(S_o) == "logical")
    S_o <- as.numeric(S_o)
  
  # Extract ML parameters
  ntree <- ml_params$ntree
  classwt_Y <- ml_params$classwt_Y
  seed <- ml_params$seed
  cores <- ml_params$cores

  # Convert R to matrix if it's a vector
  R_mat <- as.matrix(R)
  R_pred_mat <- as.matrix(R_pred)

  obs_idx <- (S_o == 1)
  if (sum(obs_idx) == 0) stop("No observations in observational sample")

  exp_idx <- (S_e == 1)
  if (sum(exp_idx) == 0) stop("No observations in experimental sample")

  # Fit E[Y | R, S_o = 1] on observations with Y available
  model_Y <- ranger::ranger(
    x = R_mat[obs_idx, , drop = FALSE],
    y = factor(Y[obs_idx], levels = c(0, 1)),
    class.weights = classwt_Y,
    num.trees = ntree,
    num.threads = cores,
    probability = TRUE,
    seed = seed
  )

  # Fit E[D | R, S_e = 1] on observations with D available
  model_D <- ranger::ranger(
    x = R_mat[exp_idx, , drop = FALSE],
    y = factor(D[exp_idx], levels = c(0, 1)),
    num.trees = ntree,
    num.threads = cores,
    probability = TRUE,
    seed = seed
  )

  # Fit P(S_e = 1 | R) on full sample
  model_S_e <- ranger::ranger(
    x = R_mat,
    y = factor(S_e, levels = c(0, 1)),
    num.trees = ntree,
    num.threads = cores,
    probability = TRUE,
    seed = seed
  )

  # Fit P(S_o = 1 | R) on full sample
  model_S_o <- ranger::ranger(
    x = R_mat,
    y = factor(S_o, levels = c(0, 1)),
    num.trees = ntree,
    num.threads = cores,
    probability = TRUE,
    seed = seed
  )

  # Initial estimate for theta
  theta_init <- get_theta_init(
    observations = data.frame(Y, D, S_e, S_o),
    predictions = data.frame(
      Y = predict(model_Y, R_mat)$predictions[,"1"],
      D = predict(model_D, R_mat)$predictions[,"1"],
      S_e = predict(model_S_e, R_mat)$predictions[,"1"],
      S_o = predict(model_S_o, R_mat)$predictions[,"1"]
    )
  )

  # Generate predictions
  predictions <- data.frame(
    Y = predict(model_Y, R_pred_mat)$predictions[,"1"],
    D = predict(model_D, R_pred_mat)$predictions[,"1"],
    S_e = predict(model_S_e, R_pred_mat)$predictions[,"1"],
    S_o = predict(model_S_o, R_pred_mat)$predictions[,"1"]
  )
  list(theta_init = theta_init, predictions = predictions)
}


#' Compute RSV estimator from predictions
#'
#' Computes the RSV treatment effect estimator given predictions and sample
#' indicators. Implements Algorithm 1 from main text with Remark 3 generalization
#' allowing observations to be in experimental only, observational only, or both.
#'
#' @param observations a data frame containing:
#'   \describe{
#'     \item{Y}{Outcome variable (binary, NA where not observed)}
#'     \item{D}{Treatment indicator (binary, NA where not observed)}
#'     \item{S_e}{Experimental sample indicator (0 or 1)}
#'     \item{S_o}{Observational sample indicator (0 or 1)}
#'   }
#' @param predictions a data frame containing:
#'   \describe{
#'     \item{Y}{Predicted \eqn{P(Y = 1 \mid R, S_o = 1)}, \eqn{\text{PRED}_Y(R)}}
#'     \item{D}{Predicted \eqn{P(D = 1 \mid R, S_e = 1)}, \eqn{\text{PRED}_D(R)}}
#'     \item{S_e}{Predicted \eqn{P(S_e = 1 \mid R)}, \eqn{\text{PRED}_{S_e}(R)}}
#'     \item{S_o}{Predicted \eqn{P(S_o = 1 \mid R)}, \eqn{\text{PRED}_{S_o}(R)}}
#'   }
#' @param theta_init Initial estimate of the treatment effect.
#' @param eps Small constant for numerical stability. Default is 1e-2.
#'
#' @return A list containing:
#'   \describe{
#'     \item{coef}{Treatment effect estimate}
#'     \item{weights}{Efficient weights used in estimation}
#'     \item{n_obs}{Sample size in observational sample}
#'     \item{n_exp}{Sample size in experimental sample}
#'     \item{n_both}{Sample size in both samples}
#'     \item{numerator}{Numerator of the treatment effect estimate}
#'     \item{denominator}{Denominator of the treatment effect estimate}
#'   }
#'
#' @keywords internal
rsv_compute <- function(observations, predictions, theta_init, eps = 1e-2) {
  if (nrow(observations) != nrow(predictions))
    stop("pred_Y, pred_D, pred_S_e, pred_S_o must have same length as Y, D, S_e, S_o")
  
  if (class(observations$S_e) == "logical")
    observations$S_e <- as.numeric(observations$S_e)
  
  if (class(observations$S_o) == "logical")
    observations$S_o <- as.numeric(observations$S_o)
  
  # Sample sizes
  n_exp <- sum(observations$S_e == 1)
  n_obs <- sum(observations$S_o == 1)
  n_both <- sum(observations$S_e == 1 & observations$S_o == 1)

  # Compute sigma^2 (Step 2d)
  sigma2 <- get_sigma2(observations = observations, predictions = predictions, theta_init = theta_init, eps = eps)

  # Efficient weight function
  Delta <- get_Delta(observations = observations, predictions = predictions)
  h <- Delta$o / sigma2

  # Compute weighted estimator
  Delta <- get_Delta_obs(observations)
  numerator <- mean(Delta$e * h)
  denominator <- mean(Delta$o * h)
  coef <- numerator / denominator

  list(
    coef = coef,
    weights = h,
    n_obs = n_obs,
    n_exp = n_exp,
    n_both = n_both,
    numerator = numerator,
    denominator = denominator,
    theta_init = theta_init
  )
}


#' Bootstrap standard errors for RSV estimator
#'
#' Computes cluster-bootstrap standard errors for the RSV estimator.
#'
#' @param observations a data frame containing:
#'   \describe{
#'     \item{Y}{Outcome variable (binary, NA where not observed)}
#'     \item{D}{Treatment indicator (binary, NA where not observed)}
#'     \item{S_e}{Experimental sample indicator (0 or 1)}
#'     \item{S_o}{Observational sample indicator (0 or 1)}
#'   }
#' @param predictions a data frame containing:
#'   \describe{
#'     \item{Y}{Predicted \eqn{P(Y = 1 \mid R, S_o = 1)}, \eqn{\text{PRED}_Y(R)}}
#'     \item{D}{Predicted \eqn{P(D = 1 \mid R, S_e = 1)}, \eqn{\text{PRED}_D(R)}}
#'     \item{S_e}{Predicted \eqn{P(S_e = 1 \mid R)}, \eqn{\text{PRED}_{S_e}(R)}}
#'     \item{S_o}{Predicted \eqn{P(S_o = 1 \mid R)}, \eqn{\text{PRED}_{S_o}(R)}}
#'   }
#' @param theta_init Initial estimate of the treatment effect.
#' @param eps Small constant for numerical stability of \code{sigma2} estimate (default \code{1e-2}).
#' @param se_params List of bootstrap parameters:
#'   \describe{
#'     \item{B}{Number of bootstrap replications (default \code{1000})}
#'     \item{fix_seed}{If \code{TRUE}, deterministic seeding is used with `set.seed(b)` for the *b*‑th replication (default \code{FALSE})}
#'     \item{cores}{Number of cores for bootstrap replications (default \code{1}).}
#'     \item{clusters}{Clusters for the bootstrap. If \code{NULL}, uses individual-level bootstrap}
#'   }
#'
#' @return A list with components:
#' \describe{
#'   \item{se}{Standard error of the treatment effect}
#'   \item{denominator_se}{Standard error of the denominator of the treatment effect}
#' }
#'
#' @keywords internal
rsv_bootstrap <- function(
    observations, predictions, theta_init, eps = 1e-2, 
    se_params = list(B = 1000, fix_seed = FALSE, cores = 1, clusters = NULL)
  ) {
  n <- length(observations$Y)
  B <- se_params$B
  fix_seed <- se_params$fix_seed
  cores <- se_params$cores
  if (is.null(se_params$clusters))
    clusters <- 1:n # If no clusters provided, use individual-level bootstrap
  else
    clusters <- se_params$clusters
  
  if (nrow(observations) != length(clusters))
    stop("pred_Y, pred_D, pred_S_e, pred_S_o must have same se_params$clusters")
  
  # Bootstrap function
  run_one <- function(b) {
    if (fix_seed)
      set.seed(b)

    # Cluster bootstrap
    unique_clusters <- unique(clusters)
    clusters_b <- sample(unique_clusters, length(unique_clusters), replace = TRUE)
    index_b <- which(clusters %in% clusters_b)
    
    observations_b <- observations[index_b, , drop = FALSE]
    predictions_b <- predictions[index_b, , drop = FALSE]
    
    out_b <- rsv_compute(
      observations = observations_b,
      predictions = predictions_b,
      theta_init = theta_init,
      eps = eps
    )
    return(c(coef = out_b$coef, denominator = out_b$denominator))
  }
  
  out <- dplyr::bind_rows(parallel::mclapply(1:B, run_one, mc.cores = cores))
  
  # Compute SE and confidence interval
  se <- sd(out$coef, na.rm = TRUE)
  denominator_se <- sd(out$denominator, na.rm = TRUE)

  list(
    se = structure(se, names = "D"),
    denominator_se = denominator_se
  )
}


#' Compute joint probabilities
#'
#' Computes joint probabilities for treatment/outcome and sample membership.
#'
#' @param x A data frame with columns \code{D}, \code{Y}, \code{S_e}, and \code{S_o}.
#'
#' @return A list containing:
#' \describe{
#'   \item{D1Se}{P(D=1, S_e=1 | R) = P(D=1 | S_e=1, R) × P(S_e=1 | R)}
#'   \item{D0Se}{P(D=0, S_e=1 | R) = P(D=0 | S_e=1, R) × P(S_e=1 | R)}
#'   \item{Y1So}{P(Y=1, S_o=1 | R) = P(Y=1 | S_o=1, R) × P(S_o=1 | R)}
#'   \item{Y0So}{P(Y=0, S_o=1 | R) = P(Y=0 | S_o=1, R) × P(S_o=1 | R)}
#' }
#'
#' @keywords internal
get_joint <- function(x) {
  list(
    D1Se =      x$D  * x$S_e, 
    D0Se = (1 - x$D) * x$S_e, 
    Y1So =      x$Y  * x$S_o, 
    Y0So = (1 - x$Y) * x$S_o  
  )
}


#' Compute marginal probability sums
#'
#' Computes sample means of joint probabilities. Used to calculate marginal
#' probabilities for treatment/outcome and sample membership.
#'
#' @param x A data frame with columns \code{D}, \code{Y}, \code{S_e}, and \code{S_o}
#'   containing observed or predicted values.
#'
#' @return A list containing mean values of:
#' \describe{
#'   \item{D1Se}{Mean of P(D=1, S_e=1 | R)}
#'   \item{D0Se}{Mean of P(D=0, S_e=1 | R)}
#'   \item{Y1So}{Mean of P(Y=1, S_o=1 | R)}
#'   \item{Y0So}{Mean of P(Y=0, S_o=1 | R)}
#' }
#'
#' @keywords internal
count_marginals <- function(observations) {
  observations$Y <- as.integer((observations$Y == 1) & (observations$S_o == 1))
  observations$D <- as.integer((observations$D == 1) & (observations$S_e == 1))
  j <- get_joint(observations)
  lapply(j, mean)
}

#' Compute treatment and outcome variations for RSV estimation (Deltas)
#' @param observations a data frame containing:
#'   \describe{
#'     \item{Y}{Outcome variable (binary, NA where not observed)}
#'     \item{D}{Treatment indicator (binary, NA where not observed)}
#'     \item{S_e}{Experimental sample indicator (0 or 1)}
#'     \item{S_o}{Observational sample indicator (0 or 1)}
#'   }
#' @param predictions a data frame containing:
#'   \describe{
#'     \item{Y}{Predicted \eqn{P(Y = 1 \mid R, S_o = 1)}, \eqn{\text{PRED}_Y(R)}}
#'     \item{D}{Predicted \eqn{P(D = 1 \mid R, S_e = 1)}, \eqn{\text{PRED}_D(R)}}
#'     \item{S_e}{Predicted \eqn{P(S_e = 1 \mid R)}, \eqn{\text{PRED}_{S_e}(R)}}
#'     \item{S_o}{Predicted \eqn{P(S_o = 1 \mid R)}, \eqn{\text{PRED}_{S_o}(R)}}
#'   }
#'
#' @return A list containing:
#' \describe{
#'   \item{e}{Treatment variation in experimental sample}
#'   \item{o}{Outcome variation in observational sample}
#' }
#'
#' @keywords internal
get_Delta <- function(observations, predictions) {
  j <- get_joint(predictions)
  count <- count_marginals(observations)
  list(
    e = j$D1Se / count$D1Se - j$D0Se / count$D0Se, # treatment variations from the Se
    o = j$Y1So / count$Y1So - j$Y0So / count$Y0So # outcome variations from the So
  )
}

#' Compute treatment and outcome variations for RSV estimation (Delta) from marginal probabilities only
#' @param observations a data frame containing:
#'   \describe{
#'     \item{Y}{Outcome variable (binary, NA where not observed)}
#'     \item{D}{Treatment indicator (binary, NA where not observed)}
#'     \item{S_e}{Experimental sample indicator (0 or 1)}
#'     \item{S_o}{Observational sample indicator (0 or 1)}
#'   }
#'
#' @return A list containing:
#' \describe{
#'   \item{e}{Treatment variation in experimental sample}
#'   \item{o}{Outcome variation in observational sample}
#' }
#'
#' @keywords internal
get_Delta_obs <- function(observations) {
  observations$Y <- as.integer((observations$Y == 1) & (observations$S_o == 1))
  observations$D <- as.integer((observations$D == 1) & (observations$S_e == 1))
  j <- get_joint(observations)
  
  count <- count_marginals(observations)
  
  list(
    e = j$D1Se / count$D1Se - j$D0Se / count$D0Se, # treatment variations from the Se
    o = j$Y1So / count$Y1So - j$Y0So / count$Y0So # outcome variations from the So
  )
}

#' Compute variance component sigma-squared
#'
#' Computes the variance component used in the efficient weight function for
#' RSV estimation. This is the conditional variance of the efficient influence
#' function.
#'
#' @param observations a data frame containing:
#'   \describe{
#'     \item{Y}{Outcome variable (binary, NA where not observed)}
#'     \item{D}{Treatment indicator (binary, NA where not observed)}
#'     \item{S_e}{Experimental sample indicator (0 or 1)}
#'     \item{S_o}{Observational sample indicator (0 or 1)}
#'   }
#' @param predictions a data frame containing:
#'   \describe{
#'     \item{Y}{Predicted \eqn{P(Y = 1 \mid R, S_o = 1)}, \eqn{\text{PRED}_Y(R)}}
#'     \item{D}{Predicted \eqn{P(D = 1 \mid R, S_e = 1)}, \eqn{\text{PRED}_D(R)}}
#'     \item{S_e}{Predicted \eqn{P(S_e = 1 \mid R)}, \eqn{\text{PRED}_{S_e}(R)}}
#'     \item{S_o}{Predicted \eqn{P(S_o = 1 \mid R)}, \eqn{\text{PRED}_{S_o}(R)}}
#'   }
#' @param theta_init Initial estimate of the treatment effect.
#' @param eps Small positive constant for numerical stability (default \code{1e-2}).
#'   Ensures sigma-squared is bounded away from zero.
#'
#' @return A numeric vector of sigma-squared values, one for each observation.
#'
#' @keywords internal
get_sigma2 <- function(observations, predictions, theta_init, eps=1e-2) {
  j <- get_joint(predictions)
  count <- count_marginals(observations)
  sigma2 <- (j$D1Se / count$D1Se^2 + j$D0Se / count$D0Se^2) + 
    theta_init^2 * (j$Y1So / count$Y1So^2 + j$Y0So / count$Y0So^2)
  
  # Lower bound on sigma for numerical stability
  sigma2 <- pmax(sigma2, eps) 
  
  return(sigma2)
}

#' Compute initial treatment effect estimate
#'
#' Computes an initial estimate of the treatment effect using predicted
#' probabilities.
#'
#' @param observations a data frame containing:
#'   \describe{
#'     \item{Y}{Outcome variable (binary, NA where not observed)}
#'     \item{D}{Treatment indicator (binary, NA where not observed)}
#'     \item{S_e}{Experimental sample indicator (0 or 1)}
#'     \item{S_o}{Observational sample indicator (0 or 1)}
#'   }
#' @param predictions a data frame containing:
#'   \describe{
#'     \item{Y}{Predicted \eqn{P(Y = 1 \mid R, S_o = 1)}, \eqn{\text{PRED}_Y(R)}}
#'     \item{D}{Predicted \eqn{P(D = 1 \mid R, S_e = 1)}, \eqn{\text{PRED}_D(R)}}
#'     \item{S_e}{Predicted \eqn{P(S_e = 1 \mid R)}, \eqn{\text{PRED}_{S_e}(R)}}
#'     \item{S_o}{Predicted \eqn{P(S_o = 1 \mid R)}, \eqn{\text{PRED}_{S_o}(R)}}
#'   }
#' @return A scalar initial treatment effect estimate.
#'
#' @keywords internal
get_theta_init <- function(observations, predictions){
  predictions$S_e <- observations$S_e
  predictions$S_o <- observations$S_o
  Delta <- get_Delta(observations, predictions)
  theta_init <- mean(Delta$e * Delta$o, na.rm = TRUE) / mean(Delta$o^2, na.rm = TRUE)
  return(theta_init)
}