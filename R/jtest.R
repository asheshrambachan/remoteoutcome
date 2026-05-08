# =============================================================================
# Hansen-Sargan J-test (overidentifying restrictions)
#
# Case (i): Incomplete observational sample, binary Y, no covariates.
#
# Two representations H_1 (efficient/RSV) and H_2 (naive) both satisfy
#   E[ H_k(R) * (Delta^e - tau * Delta^o) ] = 0
# for the same scalar tau.  Stacking 2 moments against 1 parameter gives
# one overidentifying restriction, testable as J ~ chi^2_1.
#
# Two-step GMM:
#   Step 1:  identity-weighted — tau_1 = (b'a) / (b'b)
#   Step 2:  efficient-weighted — tau_2 = (b'W a) / (b'W b),  W = Sigma^{-1}
#   J = n * g(tau_2)' W g(tau_2)   where g(tau) = a - tau * b
# =============================================================================

#' Hansen-Sargan J-test for overidentifying restrictions
#'
#' Tests whether the efficient (RSV) and naive representations yield the same
#' treatment effect estimate.  Under the identifying assumptions both moment
#' conditions hold simultaneously; the J-statistic tests this restriction.
#'
#' @param object A fitted object of class \code{"cv.rsv"}.
#' @param ... Additional arguments (unused).
#'
#' @return A list of class \code{"jtest_cv_rsv"} with elements:
#' \describe{
#'   \item{J}{J-statistic}
#'   \item{p_value}{p-value from chi-squared(1) reference distribution}
#'   \item{df}{Degrees of freedom (always 1)}
#'   \item{tau}{Two-step efficient GMM estimate of tau}
#'   \item{tau_step1}{First-step (identity-weighted) estimate of tau}
#' }
#'
#' @details
#' Requires binary Y and no covariates.  Uses cross-fitted H weights stored
#' in the \code{cv.rsv} object (both efficient and naive) together with
#' observed sample deltas Delta^e and Delta^o.
#'
#' @export
jtest <- function(object, ...) {
  UseMethod("jtest")
}

#' @export
jtest.cv.rsv <- function(object, ...) {

  y_levels <- object$y_levels
  if (length(y_levels) != 2L)
    stop("jtest() currently supports binary Y only.")
  if (!is.null(object$X_cols) && length(object$X_cols) > 0)
    stop("jtest() currently supports the no-covariate case only.")
  if (is.null(object$weights_naive))
    stop("Naive weights not found. Please refit with the current cv.rsv().")

  yK <- if (is.null(object$yK)) y_levels[1] else object$yK
  n  <- length(object$Y)

  # Build pooled df (obs + cross-fitted predictions) for computing observed deltas
  df <- cbind(
    data.frame(Y = object$Y, D = object$D, S_e = object$S_e, S_o = object$S_o),
    object$pooled_predictions
  )
  df <- .ensure_numeric_indicators(df)

  # Observed Delta^e and Delta^o  (use_pred = FALSE: uses Y, D, not predictions)
  delta_df <- compute_deltas(df, y_levels = y_levels, yK = yK,
                             X_cols = character(0), use_pred = FALSE)
  delta_e <- delta_df$delta_e
  delta_o <- as.numeric(delta_df[[grep("^delta_o_Y", names(delta_df))[1L]]])

  # H_1: efficient (RSV) weights;  H_2: naive weights  (binary Y -> single column)
  H1 <- as.numeric(object$weights[[1L]])
  H2 <- as.numeric(object$weights_naive[[1L]])

  # Drop observations with any NA in the inputs
  ok    <- complete.cases(cbind(H1, H2, delta_e, delta_o))
  H1    <- H1[ok];     H2    <- H2[ok]
  de    <- delta_e[ok]; do_   <- delta_o[ok]
  n_ok  <- sum(ok)

  if (n_ok < 3L)
    stop("Too few complete observations to compute J-statistic.")

  # g(tau) = a - tau * b  (linear in tau)
  #   a_k = (1/n) sum H_k * Delta^e
  #   b_k = (1/n) sum H_k * Delta^o
  a <- c(mean(H1 * de), mean(H2 * de))
  b <- c(mean(H1 * do_), mean(H2 * do_))

  # Step 1: identity-weighted GMM
  tau1 <- as.numeric(crossprod(b, a) / crossprod(b, b))

  # Step 2: estimate moment covariance at tau1, then efficient GMM
  u1    <- de - tau1 * do_
  M1    <- cbind(H1 * u1, H2 * u1)             # n_ok x 2
  g1    <- colMeans(M1)
  Sigma <- crossprod(sweep(M1, 2, g1)) / n_ok   # 2 x 2 sample covariance
  W     <- solve(Sigma)

  tau2  <- as.numeric((t(b) %*% W %*% a) / (t(b) %*% W %*% b))

  # J-statistic:  n * g(tau2)' W g(tau2)  ~  chi^2_1 under H0
  g2   <- a - tau2 * b
  J    <- n_ok * as.numeric(t(g2) %*% W %*% g2)
  pval <- pchisq(J, df = 1L, lower.tail = FALSE)

  out <- list(
    J          = J,
    p_value    = pval,
    df         = 1L,
    tau        = tau2,
    tau_step1  = tau1,
    n          = n_ok
  )
  class(out) <- "jtest_cv_rsv"
  out
}


#' @export
jtest.rsv <- function(object, ...) {

  y_levels <- object$y_levels
  if (length(y_levels) != 2L)
    stop("jtest() currently supports binary Y only.")
  if (!is.null(object$X_cols) && length(object$X_cols) > 0)
    stop("jtest() currently supports the no-covariate case only.")
  if (is.null(object$weights_naive))
    stop("Naive weights not found. Please refit with the current rsv() / rsv_split() / rsv_fitted().")

  yK <- if (is.null(object$yK)) y_levels[1] else object$yK
  n  <- nrow(object$data)

  df <- object$data
  df <- .ensure_numeric_indicators(df)

  # Observed Delta^e and Delta^o  (use_pred = FALSE: uses Y, D, not predictions)
  delta_df <- compute_deltas(df, y_levels = y_levels, yK = yK,
                             X_cols = character(0), use_pred = FALSE)
  delta_e <- delta_df$delta_e
  delta_o <- as.numeric(delta_df[[grep("^delta_o_Y", names(delta_df))[1L]]])

  # H_1: efficient (RSV) weights;  H_2: naive weights  (binary Y -> single column)
  H1 <- as.numeric(object$weights[[1L]])
  H2 <- as.numeric(object$weights_naive[[1L]])

  # Drop observations with any NA in the inputs
  ok    <- complete.cases(cbind(H1, H2, delta_e, delta_o))
  H1    <- H1[ok];     H2    <- H2[ok]
  de    <- delta_e[ok]; do_   <- delta_o[ok]
  n_ok  <- sum(ok)

  if (n_ok < 3L)
    stop("Too few complete observations to compute J-statistic.")

  # g(tau) = a - tau * b  (linear in tau)
  a <- c(mean(H1 * de), mean(H2 * de))
  b <- c(mean(H1 * do_), mean(H2 * do_))

  # Step 1: identity-weighted GMM
  tau1 <- as.numeric(crossprod(b, a) / crossprod(b, b))

  # Step 2: estimate moment covariance at tau1, then efficient GMM
  u1    <- de - tau1 * do_
  M1    <- cbind(H1 * u1, H2 * u1)             # n_ok x 2
  g1    <- colMeans(M1)
  Sigma <- crossprod(sweep(M1, 2, g1)) / n_ok   # 2 x 2 sample covariance
  W     <- solve(Sigma)

  tau2  <- as.numeric((t(b) %*% W %*% a) / (t(b) %*% W %*% b))

  # J-statistic:  n * g(tau2)' W g(tau2)  ~  chi^2_1 under H0
  g2   <- a - tau2 * b
  J    <- n_ok * as.numeric(t(g2) %*% W %*% g2)
  pval <- pchisq(J, df = 1L, lower.tail = FALSE)

  out <- list(
    J          = J,
    p_value    = pval,
    df         = 1L,
    tau        = tau2,
    tau_step1  = tau1,
    n          = n_ok
  )
  class(out) <- "jtest_cv_rsv"
  out
}


#' @export
print.jtest_cv_rsv <- function(x, digits = 4, ...) {
  cat("Hansen-Sargan J-test  [H_RSV vs H_naive]\n")
  cat("==========================================\n\n")
  cat(sprintf("  J-statistic: %.*f   (df = %d,  p = %.4f)\n",
              digits, x$J, x$df, x$p_value))
  cat(sprintf("  GMM tau:     %.*f\n", digits, x$tau))
  cat(sprintf("  n:           %d\n", x$n))
  invisible(x)
}
