# =============================================================================
# Relevance Extractor for RSV Objects
#
# relevance() extracts E_n[ H(R) Delta^o ] from a fitted cv.rsv object,
# combining both the RSV and naive estimators alongside bootstrap SEs and
# normal CIs.
# =============================================================================

#' Extract relevance table from a fitted RSV object
#'
#' Returns a data frame with the relevance \eqn{E_n[H_d(R) \cdot \Delta^o(d)]} for
#' both the RSV and naive estimators, along with bootstrap SEs and normal CIs.
#'
#' @param object A fitted object of class \code{"cv.rsv"}.
#' @param level Confidence level for normal CI (default 0.95).
#' @param ... Additional arguments (unused).
#'
#' @return A data frame of class \code{"relevance_table"} with columns:
#' \describe{
#'   \item{estimator}{"rsv" or "naive"}
#'   \item{estimate}{Point estimate of \eqn{E_n[H \cdot \delta_o]}}
#'   \item{se}{Bootstrap SE (NA if \code{add_se} has not been run)}
#'   \item{ci_lower}{Lower normal CI bound}
#'   \item{ci_upper}{Upper normal CI bound}
#' }
#'
#' @export
relevance <- function(object, level = 0.95, ...) {
  UseMethod("relevance")
}

#' @export
relevance.rsv <- function(object, level = 0.95, ...) {
  if (is.null(object$relevance))
    stop("Relevance estimates missing. Please refit with the current rsv() / rsv_split() / rsv_fitted().")

  z <- qnorm(1 - (1 - level) / 2)

  make_row <- function(estimator, est, se) {
    data.frame(
      estimator = estimator,
      estimate  = est,
      se        = se,
      ci_lower  = if (is.na(se)) NA_real_ else est - z * se,
      ci_upper  = if (is.na(se)) NA_real_ else est + z * se,
      stringsAsFactors = FALSE
    )
  }

  rows <- list(
    make_row("rsv",   object$relevance,       object$relevance_se   %||% NA_real_),
    make_row("naive", object$relevance_naive,  object$relevance_naive_se %||% NA_real_)
  )

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "level") <- level
  attr(out, "B")     <- object$bootstrap$B %||% NA
  class(out) <- c("relevance_table", "data.frame")
  out
}


#' @export
relevance.cv.rsv <- function(object, level = 0.95, ...) {
  if (is.null(object$relevance))
    stop("Relevance estimates missing. Please refit with the current cv.rsv().")

  z <- qnorm(1 - (1 - level) / 2)

  make_row <- function(estimator, est, se) {
    data.frame(
      estimator = estimator,
      estimate  = est,
      se        = se,
      ci_lower  = if (is.na(se)) NA_real_ else est - z * se,
      ci_upper  = if (is.na(se)) NA_real_ else est + z * se,
      stringsAsFactors = FALSE
    )
  }

  rows <- list(
    make_row("rsv",   object$relevance,       object$relevance_se   %||% NA_real_),
    make_row("naive", object$relevance_naive,  object$relevance_naive_se %||% NA_real_)
  )

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "level") <- level
  attr(out, "B")     <- object$bootstrap$B %||% NA
  class(out) <- c("relevance_table", "data.frame")
  out
}


#' @export
print.relevance_table <- function(x, digits = 4, ...) {
  level <- attr(x, "level") %||% 0.95
  B     <- attr(x, "B")     %||% NA
  if (!is.na(B)) {
    cat(sprintf("Relevance  E_n[ H(R) Delta^o ]   (B=%d, %.0f%% normal CI)\n",
                B, 100 * level))
  } else {
    cat(sprintf("Relevance  E_n[ H(R) Delta^o ]   (%.0f%% normal CI)\n",
                100 * level))
  }
  cat(strrep("-", 60), "\n")
  fmt           <- x
  num_cols      <- c("estimate", "se", "ci_lower", "ci_upper")
  fmt[num_cols] <- lapply(fmt[num_cols], round, digits = digits)
  print.data.frame(fmt, row.names = FALSE)
  invisible(x)
}
