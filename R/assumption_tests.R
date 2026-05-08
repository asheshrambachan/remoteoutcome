# =============================================================================
# Assumption Tests for the RSV Estimator
# =============================================================================

#' Test the stability assumption
#'
#' Tests Assumption 2(i): \eqn{R \perp S \mid Y}, i.e. the sensing mechanism
#' \eqn{f(R \mid Y)} is the same in the experimental and observational samples.
#' For each outcome level \eqn{y} and each RSV column, compares the distribution
#' of R in \eqn{\{S_e = 1, Y = y\}} against \eqn{\{S_o = 1, Y = y\}} via a
#' two-sample Kolmogorov-Smirnov test.
#'
#' Failure to reject supports stability: the R--Y relationship learned in
#' one sample can be transported to the other.
#'
#' @param R Remote sensing variables (numeric matrix or data frame, \code{n x p}).
#' @param Y Outcome vector. Must be observed (non-NA) in at least one of the
#'   two samples for the test to be informative.
#' @param S_e Experimental sample indicator (0/1 or logical, length \code{n}).
#' @param S_o Observational sample indicator (0/1 or logical, length \code{n}).
#' @param y_levels Levels of Y to test (default: all unique non-NA values, sorted).
#'
#' @return A data frame of class \code{"stability_test"} with one row per
#'   (Y level, R column) combination and columns:
#' \describe{
#'   \item{y_level}{Outcome level being tested}
#'   \item{r_col}{RSV column name}
#'   \item{n_exp}{Number of experimental units at this Y level}
#'   \item{n_obs}{Number of observational units at this Y level}
#'   \item{ks_stat}{KS test statistic}
#'   \item{p_value}{Two-sample KS p-value}
#' }
#' Rows where either group has fewer than 2 observations have \code{NA}
#' statistics.
#'
#' @details
#' This test requires Y to be observed in **both** samples. It cannot be run
#' in a pure "no experimental outcomes" design where Y is only measured in the
#' observational sample — which is exactly the main use case for the RSV
#' estimator. In that setting, use the J-test (\code{\link{jtest}}) as an
#' indirect check.
#'
#' A KS test is a necessary but not sufficient check: it tests the marginal
#' distribution of each R column conditional on Y, not the full joint
#' conditional. Low power is expected at small sample sizes.
#'
#' @seealso \code{\link{test_no_direct_effect}}, \code{\link{jtest}}
#'
#' @examples
#' # n_v = 100 overlap units (S_e = 1, S_o = 1) have Y observed
#' dat <- sim_rsv_data(n_e = 300, n_o = 700, n_v = 100, seed = 42)
#' R   <- as.matrix(dat[, paste0("R", 1:5)])
#' res <- test_stability(R, dat$Y, dat$S_e, dat$S_o)
#' print(res)
#'
#' @export
test_stability <- function(R, Y, S_e, S_o, y_levels = NULL) {
  cl <- match.call()
  n  <- length(Y)

  R   <- as.matrix(R)
  S_e <- as.integer(S_e)
  S_o <- as.integer(S_o)

  if (nrow(R) != n)  stop("R must have the same number of rows as Y.")
  if (length(S_e) != n) stop("S_e must have the same length as Y.")
  if (length(S_o) != n) stop("S_o must have the same length as Y.")

  if (all(is.na(Y[S_e == 1L])))
    stop(
      "Y is not observed in the experimental sample (S_e = 1).\n",
      "The stability assumption cannot be tested without experimental outcomes.\n",
      "Use jtest() as an indirect check instead."
    )
  if (all(is.na(Y[S_o == 1L])))
    stop("Y is not observed in the observational sample (S_o = 1).")

  if (is.null(y_levels)) y_levels <- sort(unique(na.omit(Y)))

  r_cols <- if (!is.null(colnames(R))) colnames(R) else paste0("R", seq_len(ncol(R)))
  colnames(R) <- r_cols

  rows <- vector("list", length(y_levels) * ncol(R))
  idx  <- 0L
  for (y_val in y_levels) {
    mask_exp <- S_e == 1L & !is.na(Y) & Y == y_val
    mask_obs <- S_o == 1L & !is.na(Y) & Y == y_val
    for (j in seq_len(ncol(R))) {
      idx <- idx + 1L
      r_e <- R[mask_exp, j]
      r_o <- R[mask_obs, j]
      if (length(r_e) < 2L || length(r_o) < 2L) {
        rows[[idx]] <- data.frame(
          y_level = y_val, r_col = r_cols[j],
          n_exp = length(r_e), n_obs = length(r_o),
          ks_stat = NA_real_, p_value = NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        ks <- stats::ks.test(r_e, r_o)
        rows[[idx]] <- data.frame(
          y_level = y_val, r_col = r_cols[j],
          n_exp = length(r_e), n_obs = length(r_o),
          ks_stat = unname(ks$statistic),
          p_value = ks$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "call") <- cl
  class(out) <- c("stability_test", "data.frame")
  out
}


#' Test the no-direct-effect assumption
#'
#' Tests Assumption 3(ii): \eqn{R \perp D \mid S_e = 1, Y}, i.e. treatment
#' affects Y but does not independently shift R conditional on Y.
#' For each outcome level \eqn{y} and each RSV column, compares the
#' distribution of R among treated (\eqn{D = 1}) vs untreated (\eqn{D = 0})
#' experimental units with \eqn{Y = y} via a two-sample KS test.
#'
#' Failure to reject supports the no-direct-effect assumption: R can serve
#' as a valid proxy for Y in the observational sample.
#'
#' @param R Remote sensing variables (numeric matrix or data frame, \code{n x p}).
#' @param Y Outcome vector. Must be observed (non-NA) in the experimental
#'   sample for the test to be run.
#' @param D Treatment indicator (0/1 or logical, length \code{n}).
#'   \code{NA} is allowed outside the experimental sample.
#' @param S_e Experimental sample indicator (0/1 or logical, length \code{n}).
#' @param y_levels Levels of Y to test (default: all unique non-NA values, sorted).
#'
#' @return A data frame of class \code{"no_direct_effect_test"} with one row
#'   per (Y level, R column) combination and columns:
#' \describe{
#'   \item{y_level}{Outcome level being tested}
#'   \item{r_col}{RSV column name}
#'   \item{n_d0}{Number of untreated experimental units at this Y level}
#'   \item{n_d1}{Number of treated experimental units at this Y level}
#'   \item{ks_stat}{KS test statistic}
#'   \item{p_value}{Two-sample KS p-value}
#' }
#'
#' @details
#' This test requires Y to be observed in the experimental sample **and** both
#' treatment arms to be present. It cannot be run in the typical RSV setting
#' where Y is not measured for experimental units (the main use case of the
#' estimator). Use \code{\link{jtest}} as an indirect check in that case.
#'
#' @seealso \code{\link{test_stability}}, \code{\link{jtest}}
#'
#' @examples
#' # n_v = 100 overlap units (S_e = 1, S_o = 1) have Y and D observed
#' dat <- sim_rsv_data(n_e = 300, n_o = 700, n_v = 100, seed = 42)
#' R   <- as.matrix(dat[, paste0("R", 1:5)])
#' res <- test_no_direct_effect(R, dat$Y, dat$D, dat$S_e)
#' print(res)
#'
#' @export
test_no_direct_effect <- function(R, Y, D, S_e, y_levels = NULL) {
  cl <- match.call()
  n  <- length(Y)

  R   <- as.matrix(R)
  S_e <- as.integer(S_e)
  D   <- as.integer(D)

  if (nrow(R) != n)     stop("R must have the same number of rows as Y.")
  if (length(S_e) != n) stop("S_e must have the same length as Y.")
  if (length(D) != n)   stop("D must have the same length as Y.")

  exp_mask <- S_e == 1L

  if (all(is.na(Y[exp_mask])))
    stop(
      "Y is not observed in the experimental sample (S_e = 1).\n",
      "The no-direct-effect assumption cannot be tested without experimental outcomes.\n",
      "Use jtest() as an indirect check instead."
    )

  d_exp <- D[exp_mask]
  if (!any(d_exp == 0L, na.rm = TRUE) || !any(d_exp == 1L, na.rm = TRUE))
    stop("Both D = 0 and D = 1 must be present in the experimental sample.")

  if (is.null(y_levels)) y_levels <- sort(unique(na.omit(Y)))

  r_cols <- if (!is.null(colnames(R))) colnames(R) else paste0("R", seq_len(ncol(R)))
  colnames(R) <- r_cols

  rows <- vector("list", length(y_levels) * ncol(R))
  idx  <- 0L
  for (y_val in y_levels) {
    mask_d0 <- exp_mask & !is.na(D) & D == 0L & !is.na(Y) & Y == y_val
    mask_d1 <- exp_mask & !is.na(D) & D == 1L & !is.na(Y) & Y == y_val
    for (j in seq_len(ncol(R))) {
      idx <- idx + 1L
      r_d0 <- R[mask_d0, j]
      r_d1 <- R[mask_d1, j]
      if (length(r_d0) < 2L || length(r_d1) < 2L) {
        rows[[idx]] <- data.frame(
          y_level = y_val, r_col = r_cols[j],
          n_d0 = length(r_d0), n_d1 = length(r_d1),
          ks_stat = NA_real_, p_value = NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        ks <- stats::ks.test(r_d0, r_d1)
        rows[[idx]] <- data.frame(
          y_level = y_val, r_col = r_cols[j],
          n_d0 = length(r_d0), n_d1 = length(r_d1),
          ks_stat = unname(ks$statistic),
          p_value = ks$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "call") <- cl
  class(out) <- c("no_direct_effect_test", "data.frame")
  out
}


#' @export
print.stability_test <- function(x, digits = 3, ...) {
  cat("Stability assumption test  [F(R | S=e, Y=y) = F(R | S=o, Y=y)]\n")
  cat(strrep("-", 62), "\n")
  .print_ks_table(x, digits)
  cat("\nH0: distributions are equal (KS test).\n")
  cat("Failure to reject supports stability.\n")
  invisible(x)
}


#' @export
print.no_direct_effect_test <- function(x, digits = 3, ...) {
  cat("No-direct-effect assumption test  [F(R | S=e, D=0, Y=y) = F(R | S=e, D=1, Y=y)]\n")
  cat(strrep("-", 70), "\n")
  .print_ks_table(x, digits)
  cat("\nH0: distributions are equal (KS test).\n")
  cat("Failure to reject supports no direct effect of treatment on R.\n")
  invisible(x)
}


# Shared table printer for both test classes
.print_ks_table <- function(x, digits) {
  fmt           <- x
  fmt$ks_stat   <- round(fmt$ks_stat, digits)
  fmt$p_value   <- round(fmt$p_value, digits)
  print.data.frame(fmt, row.names = FALSE)
}
