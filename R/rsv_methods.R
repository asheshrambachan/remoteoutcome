# =============================================================================
# S3 Methods for rsv class
# =============================================================================

#' Print method for rsv objects
#'
#' @param x An object of class "rsv".
#' @param ... Additional arguments (unused).
#'
#' @export
print.rsv <- function(x, ...) {
  cat("RSV Treatment Effect Estimate\n")
  cat("==============================\n\n")
  cat(sprintf("Coefficient: %.4f", x$coef))
  if (!is.null(x$se)) {
    cat(sprintf(" (SE: %.4f)\n", x$se))
    ci_level <- (1 - x$alpha) * 100
    cat(sprintf("%.0f%% CI: [%.4f, %.4f]\n", ci_level, x$ci_lower, x$ci_upper))
  } else {
    cat("\n")
  }
  cat("\n")
  cat(sprintf("Sample sizes:\n"))
  cat(sprintf("  Experimental: %d\n", x$n_exp))
  cat(sprintf("  Observational: %d\n", x$n_obs))
  if (!is.null(x$n_both) && x$n_both > 0) {
    cat(sprintf("  Both: %d\n", x$n_both))
  }
  if (!is.null(x$method)) {
    cat(sprintf("\nMethod: %s\n", x$method))
  }
  invisible(x)
}


#' Summary method for rsv objects
#'
#' @param object An object of class "rsv".
#' @param ... Additional arguments (unused).
#'
#' @export
summary.rsv <- function(object, ...) {
  cat("RSV Treatment Effect Estimate\n")
  cat("==============================\n\n")

  # Create coefficient table
  if (!is.null(object$se)) {
    t_stat <- object$coef / object$se
    p_value <- 2 * (1 - pnorm(abs(t_stat)))

    coef_table <- data.frame(
      Estimate = object$coef,
      Std.Error = object$se,
      t.value = t_stat,
      Pr..t. = p_value,
      CI.Lower = object$ci_lower,
      CI.Upper = object$ci_upper
    )
    rownames(coef_table) <- "Treatment Effect"

    cat("Coefficient:\n")
    print(coef_table, digits = 4)
    cat("\n")
  } else {
    cat(sprintf("Coefficient: %.4f\n\n", object$coef))
  }

  # Sample information
  cat("Sample sizes:\n")
  cat(sprintf("  Experimental: %d\n", object$n_exp))
  cat(sprintf("  Observational: %d\n", object$n_obs))
  if (!is.null(object$n_both) && object$n_both > 0) {
    cat(sprintf("  Both: %d\n", object$n_both))
  }
  cat("\n")

  # Method information
  if (!is.null(object$method)) {
    cat(sprintf("Prediction fitting method: %s\n", object$method))
  }

  # Call information
  if (!is.null(object$call)) {
    cat("\nCall:\n")
    print(object$call)
  }

  invisible(object)
}


#' Extract coefficients from rsv objects
#'
#' @param object An object of class "rsv".
#' @param ... Additional arguments (unused).
#'
#' @return The treatment effect coefficient.
#'
#' @export
coef.rsv <- function(object, ...) {
  structure(
    object$coef,
    names = "Treatment Effect"
  )
}


#' Extract variance-covariance matrix from rsv objects
#'
#' @param object An object of class "rsv".
#' @param ... Additional arguments (unused).
#'
#' @return A 1x1 matrix containing the variance of the treatment effect.
#'
#' @export
vcov.rsv <- function(object, ...) {
  if (is.null(object$se)) {
    stop("Standard error not available. Re-run with se = TRUE.")
  }
  structure(
    matrix(object$se^2, nrow = 1, ncol = 1),
    dimnames = list("Treatment Effect", "Treatment Effect")
  )
}


#' Confidence intervals for rsv objects
#'
#' @param object An object of class "rsv".
#' @param parm Parameter (ignored, included for S3 compatibility).
#' @param level Confidence level (default 0.95). Note: this returns the CI
#'   computed during estimation based on the alpha parameter. To change the
#'   confidence level, re-run \code{rsv_estimate()} with a different alpha.
#' @param ... Additional arguments (unused).
#'
#' @return A matrix with lower and upper confidence bounds.
#'
#' @export
confint.rsv <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$ci_lower) || is.null(object$ci_upper)) {
    stop("Confidence intervals not available. Re-run with se = TRUE.")
  }

  # Check if requested level matches computed level
  computed_level <- 1 - object$alpha
  if (abs(level - computed_level) > 1e-6) {
    warning(sprintf(
      "Returning %.0f%% CI computed during estimation. To get %.0f%% CI, re-run rsv_estimate() with alpha = %.4f.",
      computed_level * 100, level * 100, 1 - level
    ))
  }

  # Determine column names based on alpha
  lower_pct <- sprintf("%.1f %%", object$alpha / 2 * 100)
  upper_pct <- sprintf("%.1f %%", (1 - object$alpha / 2) * 100)

  structure(
    matrix(c(object$ci_lower, object$ci_upper), nrow = 1, ncol = 2),
    dimnames = list("Treatment Effect", c(lower_pct, upper_pct))
  )
}
