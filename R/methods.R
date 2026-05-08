# =============================================================================
# S3 Methods for RSV Objects
# =============================================================================

# =============================================================================
# Methods for class "rsv"
# =============================================================================

#' Print method for rsv objects
#'
#' @param x An object of class "rsv"
#' @param digits Number of digits to print (default 4)
#' @param ... Additional arguments (unused)
#'
#' @export
print.rsv <- function(x, digits = 4, ...) {
  cat("RSV Treatment Effect Estimate\n")
  cat("==============================\n\n")

  cat("RSV estimator\n")
  cat(sprintf("  ATE: %.*f", digits, x$coefficients))
  if (!is.null(x$se)) cat(sprintf(" (SE: %.*f)", digits, x$se))
  cat("\n")
  if (!is.null(x$relevance) && !is.na(x$relevance)) cat(sprintf("  Relevance: %.*f", digits, x$relevance))
  if (!is.null(x$relevance_se)) cat(sprintf(" (SE: %.*f)", digits, x$relevance_se))
  if (!is.null(x$relevance) && !is.na(x$relevance)) cat("\n")

  if (!is.null(x$coefficients_naive)) {
    cat("\nNaive estimator  [H = P(Y=1|R)]\n")
    cat(sprintf("  ATE: %.*f", digits, x$coefficients_naive))
    if (!is.null(x$se_naive)) cat(sprintf(" (SE: %.*f)", digits, x$se_naive))
    cat("\n")
    if (!is.null(x$relevance_naive)) cat(sprintf("  Relevance: %.*f", digits, x$relevance_naive))
    if (!is.null(x$relevance_naive_se)) cat(sprintf(" (SE: %.*f)", digits, x$relevance_naive_se))
    if (!is.null(x$relevance_naive)) cat("\n")
  }

  cat(sprintf("\nSample sizes:\n"))
  cat(sprintf("  Experimental: %d\n", x$n.exp))
  cat(sprintf("  Observational: %d\n", x$n.obs))
  if (!is.null(x$n.both) && x$n.both > 0) {
    cat(sprintf("  Both: %d\n", x$n.both))
  }

  cat(sprintf("\nMethod: %s", x$split))
  if (x$split == "sample" && !is.null(x$train.proportion)) {
    cat(sprintf(" (train: %.0f%%, test: %d obs)",
                x$train.proportion * 100, x$test.size))
  } else if (x$split == "fitted" && !is.null(x$test.size)) {
    cat(sprintf(" (%d obs used)", x$test.size))
  }
  cat("\n")

  # Display model types
  if (!is.null(x$models)) {
    model_types <- sapply(x$models, function(m) m$model)
    if (length(unique(model_types)) == 1) {
      cat(sprintf("Models: all %s\n", model_types[1]))
    } else {
      cat("Models:\n")
      for (nm in names(model_types)) {
        cat(sprintf("  %s: %s\n", nm, model_types[nm]))
      }
    }
  }

  invisible(x)
}


#' Summary method for rsv objects
#'
#' @param object An object of class "rsv"
#' @param ... Additional arguments (unused)
#'
#' @export
summary.rsv <- function(object, ...) {
  cat("RSV Treatment Effect Estimate\n")
  cat("==============================\n\n")

  # Coefficient table
  if (!is.null(object$se)) {
    t_stat <- object$coefficients / object$se
    p_value <- 2 * (1 - stats::pnorm(abs(t_stat)))

    coef_table <- cbind(
      Estimate = object$coefficients,
      "Std. Error" = object$se,
      "z value" = t_stat,
      "Pr(>|z|)" = p_value
    )
    rownames(coef_table) <- "D"

    cat("Coefficients:\n")
    stats::printCoefmat(coef_table, digits = 4, signif.stars = TRUE,
                        P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  } else {
    cat(sprintf("Coefficient: %.4f\n", object$coefficients))
    cat("(Standard errors not computed. Use add_se() to add them.)\n\n")
  }

  # Sample information
  cat("Sample sizes:\n")
  cat(sprintf("  Experimental: %d\n", object$n.exp))
  cat(sprintf("  Observational: %d\n", object$n.obs))
  if (!is.null(object$n.both) && object$n.both > 0) {
    cat(sprintf("  Both: %d\n", object$n.both))
  }
  cat("\n")

  # Method information
  cat(sprintf("Method: %s\n", object$split))

  # Display model types
  if (!is.null(object$models)) {
    model_types <- sapply(object$models, function(m) m$model)
    if (length(unique(model_types)) == 1) {
      cat(sprintf("Models: all %s\n", model_types[1]))
    } else {
      cat("Models:\n")
      for (nm in names(model_types)) {
        cat(sprintf("  %s: %s\n", nm, model_types[nm]))
      }
    }
  }

  # Call information
  if (!is.null(object$call)) {
    cat("\nCall:\n")
    print(object$call)
  }

  if (!is.null(object$se.call)) {
    cat("\nSE Call:\n")
    print(object$se.call)
  }

  invisible(object)
}


#' Extract coefficients from rsv objects
#'
#' @param object An object of class "rsv"
#' @param ... Additional arguments (unused)
#'
#' @return Named numeric vector of coefficients
#'
#' @export
coef.rsv <- function(object, ...) {
  object$coefficients
}


#' Extract variance-covariance matrix from rsv objects
#'
#' @param object An object of class "rsv"
#' @param ... Additional arguments (unused)
#'
#' @return 1x1 variance-covariance matrix
#'
#' @export
vcov.rsv <- function(object, ...) {
  if (is.null(object$se)) {
    stop("Standard errors not available. Use add_se() to compute them.")
  }

  if (!is.null(object$vcov)) {
    return(object$vcov)
  }

  # Construct vcov from se
  structure(
    matrix(object$se^2, 1, 1),
    dimnames = list("D", "D")
  )
}


#' Confidence intervals for rsv objects
#'
#' @param object An object of class "rsv"
#' @param parm Parameter names (default "D")
#' @param level Confidence level (default 0.95)
#' @param ... Additional arguments (unused)
#'
#' @return Matrix with confidence intervals
#'
#' @export
confint.rsv <- function(object, parm = "D", level = 0.95, ...) {

  if (is.null(object$se)) {
    stop("Standard errors not available. Use add_se() to compute them.")
  }

  if (parm != "D") {
    stop("Only parameter 'D' is available for rsv objects")
  }

  alpha <- 1 - level
  z_crit <- stats::qnorm(1 - alpha / 2)

  ci_lower <- object$coefficients - z_crit * object$se
  ci_upper <- object$coefficients + z_crit * object$se

  lower_pct <- sprintf("%.1f %%", alpha / 2 * 100)
  upper_pct <- sprintf("%.1f %%", (1 - alpha / 2) * 100)

  structure(
    matrix(c(ci_lower, ci_upper), nrow = 1, ncol = 2),
    dimnames = list(parm, c(lower_pct, upper_pct))
  )
}


# =============================================================================
# Methods for class "cv.rsv"
# =============================================================================

#' Print method for cv.rsv objects
#'
#' @param x An object of class "cv.rsv"
#' @param digits Number of digits to print (default 4)
#' @param ... Additional arguments (unused)
#'
#' @export
print.cv.rsv <- function(x, digits = 4, ...) {
  cat("RSV Treatment Effect Estimate (Cross-Fitted)\n")
  cat("=============================================\n\n")

  cat("RSV estimator\n")
  cat(sprintf("  ATE: %.*f", digits, x$coefficients))
  if (!is.null(x$se)) cat(sprintf(" (SE: %.*f)", digits, x$se))
  cat("\n")
  if (!is.null(x$cond_J))
    cat(sprintf("  Score matrix kappa(J): %.1f\n", x$cond_J))
  if (!is.null(x$cond_J_median) && !is.na(x$cond_J_median))
    cat(sprintf("  Score matrix kappa(J): median=%.1f  max=%.1f  (over %d bootstrap draws)\n",
                x$cond_J_median, x$cond_J_max, x$bootstrap$B))
  if (!is.null(x$relevance)) cat(sprintf("  Relevance: %.*f", digits, x$relevance))
  if (!is.null(x$relevance_se)) cat(sprintf(" (SE: %.*f)", digits, x$relevance_se))
  if (!is.null(x$relevance)) cat("\n")

  if (!is.null(x$coefficients_naive)) {
    cat("\nNaive estimator  [H = P(Y=1|R)]\n")
    cat(sprintf("  ATE: %.*f", digits, x$coefficients_naive))
    if (!is.null(x$se_naive)) cat(sprintf(" (SE: %.*f)", digits, x$se_naive))
    cat("\n")
    if (!is.null(x$relevance_naive)) cat(sprintf("  Relevance: %.*f", digits, x$relevance_naive))
    if (!is.null(x$relevance_naive_se)) cat(sprintf(" (SE: %.*f)", digits, x$relevance_naive_se))
    if (!is.null(x$relevance_naive)) cat("\n")
  }

  cat(sprintf("\nSample sizes:\n"))
  cat(sprintf("  Experimental: %d\n", x$n.exp))
  cat(sprintf("  Observational: %d\n", x$n.obs))
  if (!is.null(x$n.both) && x$n.both > 0) {
    cat(sprintf("  Both: %d\n", x$n.both))
  }

  cat(sprintf("\nMethod: %d-fold cross-fitting\n", x$nfolds))

  # Display model types
  if (!is.null(x$models)) {
    model_types <- sapply(x$models, function(m) m$model)
    if (length(unique(model_types)) == 1) {
      cat(sprintf("Models: all %s\n", model_types[1]))
    } else {
      cat("Models:\n")
      for (nm in names(model_types)) {
        cat(sprintf("  %s: %s\n", nm, model_types[nm]))
      }
    }
  }

  invisible(x)
}


#' Summary method for cv.rsv objects
#'
#' @param object An object of class "cv.rsv"
#' @param ... Additional arguments (unused)
#'
#' @export
summary.cv.rsv <- function(object, ...) {
  cat("RSV Treatment Effect Estimate (Cross-Fitted)\n")
  cat("=============================================\n\n")

  # Coefficient table
  if (!is.null(object$se)) {
    t_stat <- object$coefficients / object$se
    p_value <- 2 * (1 - stats::pnorm(abs(t_stat)))

    coef_table <- cbind(
      Estimate = object$coefficients,
      "Std. Error" = object$se,
      "z value" = t_stat,
      "Pr(>|z|)" = p_value
    )
    rownames(coef_table) <- "D"

    cat("Average coefficient:\n")
    stats::printCoefmat(coef_table, digits = 4, signif.stars = TRUE,
                        P.values = TRUE, has.Pvalue = TRUE)
    if (!is.null(object$cond_J))
      cat(sprintf("Score matrix kappa(J): %.1f\n", object$cond_J))
    if (!is.null(object$cond_J_median) && !is.na(object$cond_J_median))
      cat(sprintf("Score matrix kappa(J): median=%.1f  max=%.1f  (over %d bootstrap draws)\n",
                  object$cond_J_median, object$cond_J_max, object$bootstrap$B))
    cat("\n")
  } else {
    cat(sprintf("Average coefficient: %.4f\n", object$coefficients))
    cat("(Standard errors not computed. Use add_se() to add them.)\n\n")
  }

  cat("\n")

  # Sample information
  cat("Sample sizes:\n")
  cat(sprintf("  Experimental: %d\n", object$n.exp))
  cat(sprintf("  Observational: %d\n", object$n.obs))
  if (!is.null(object$n.both) && object$n.both > 0) {
    cat(sprintf("  Both: %d\n", object$n.both))
  }
  cat("\n")

  cat(sprintf("Method: %d-fold cross-fitting\n", object$nfolds))

  # Display model types
  if (!is.null(object$models)) {
    model_types <- sapply(object$models, function(m) m$model)
    if (length(unique(model_types)) == 1) {
      cat(sprintf("Models: all %s\n", model_types[1]))
    } else {
      cat("Models:\n")
      for (nm in names(model_types)) {
        cat(sprintf("  %s: %s\n", nm, model_types[nm]))
      }
    }
  }

  # Call information
  if (!is.null(object$call)) {
    cat("\nCall:\n")
    print(object$call)
  }

  if (!is.null(object$se.call)) {
    cat("\nSE Call:\n")
    print(object$se.call)
  }

  invisible(object)
}


#' Extract coefficients from cv.rsv objects
#'
#' @param object An object of class "cv.rsv"
#' @param ... Additional arguments (unused)
#'
#' @return Named numeric vector of average coefficient
#'
#' @export
coef.cv.rsv <- function(object, ...) {
  object$coefficients
}


#' Extract variance-covariance matrix from cv.rsv objects
#'
#' @param object An object of class "cv.rsv"
#' @param ... Additional arguments (unused)
#'
#' @return 1x1 variance-covariance matrix
#'
#' @export
vcov.cv.rsv <- function(object, ...) {
  if (is.null(object$se)) {
    stop("Standard errors not available. Use add_se() to compute them.")
  }

  if (!is.null(object$vcov)) {
    return(object$vcov)
  }

  structure(
    matrix(object$se^2, 1, 1),
    dimnames = list("D", "D")
  )
}


#' Confidence intervals for cv.rsv objects
#'
#' @param object An object of class "cv.rsv"
#' @param parm Parameter names (default "D")
#' @param level Confidence level (default 0.95)
#' @param ... Additional arguments (unused)
#'
#' @return Matrix with confidence intervals
#'
#' @export
confint.cv.rsv <- function(object, parm = "D", level = 0.95, ...) {

  if (is.null(object$se)) {
    stop("Standard errors not available. Use add_se() to compute them.")
  }

  if (parm != "D") {
    stop("Only parameter 'D' is available for cv.rsv objects")
  }

  alpha <- 1 - level
  z_crit <- stats::qnorm(1 - alpha / 2)

  ci_lower <- object$coefficients - z_crit * object$se
  ci_upper <- object$coefficients + z_crit * object$se

  lower_pct <- sprintf("%.1f %%", alpha / 2 * 100)
  upper_pct <- sprintf("%.1f %%", (1 - alpha / 2) * 100)

  structure(
    matrix(c(ci_lower, ci_upper), nrow = 1, ncol = 2),
    dimnames = list(parm, c(lower_pct, upper_pct))
  )
}


