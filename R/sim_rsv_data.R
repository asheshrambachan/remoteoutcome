# =============================================================================
# Simulate RSV Data
# =============================================================================

#' Simulate a dataset for the RSV estimator
#'
#' Generates a synthetic dataset with an experimental sample, an optional
#' overlap (validation) subsample, and an observational sample, suitable for
#' illustrating and testing \code{\link{cv.rsv}} and related functions.
#' The outcome is binary (K = 2).
#'
#' @param n_e Number of experimental units (default 300). Includes any
#'   overlap units (\code{n_v}).
#' @param n_o Number of purely observational units (default 700).
#' @param n_v Number of overlap units (default 0): units in both the
#'   experimental and observational samples (\code{S_e = 1, S_o = 1}).
#'   Must satisfy \code{0 <= n_v <= n_e}. Overlap units have both \code{D}
#'   and \code{Y} observed; the remaining \code{n_e - n_v} purely experimental
#'   units have \code{Y = NA}.
#' @param tau True average treatment effect (default 0.10).
#' @param p0 Baseline P(Y = 1) in the observational sample (default 0.30).
#' @param n_r Number of remotely sensed variable (RSV) columns (default 5).
#' @param seed Random seed for reproducibility (default \code{NULL}).
#'
#' @return A data frame with \code{n_e + n_o} rows and columns:
#' \describe{
#'   \item{Y}{Binary outcome (0/1). \code{NA} for purely experimental units
#'     (\code{S_e = 1, S_o = 0}).}
#'   \item{D}{Binary treatment (0/1). \code{NA} for purely observational units
#'     (\code{S_e = 0, S_o = 1}).}
#'   \item{S_e}{Experimental sample indicator (1 = in experimental sample).}
#'   \item{S_o}{Observational sample indicator (1 = in observational sample).}
#'   \item{R1, R2, ...}{Remotely sensed variables, correlated with true Y.}
#' }
#'
#' @details
#' The dataset has three types of units:
#' \itemize{
#'   \item Pure experimental (\code{S_e = 1, S_o = 0}, \code{n_e - n_v} units):
#'     \eqn{D \sim \mathrm{Bernoulli}(0.5)}, \eqn{Y = \mathrm{NA}}.
#'   \item Overlap (\code{S_e = 1, S_o = 1}, \code{n_v} units):
#'     \eqn{D \sim \mathrm{Bernoulli}(0.5)},
#'     \eqn{Y | D=0 \sim \mathrm{Bernoulli}(p_0)},
#'     \eqn{Y | D=1 \sim \mathrm{Bernoulli}(p_0 + \tau)}.
#'   \item Pure observational (\code{S_e = 0, S_o = 1}, \code{n_o} units):
#'     \eqn{Y \sim \mathrm{Bernoulli}(p_0)}, \eqn{D = \mathrm{NA}}.
#' }
#' RSVs are generated as
#' \eqn{R_k \sim \mathcal{N}(0.5 \cdot Y_{\mathrm{true}}, 1)} for each column,
#' using the latent Y for pure experimental units.
#'
#' When \code{n_v = 0} (the default), Y is entirely unobserved in the
#' experimental sample — the main use case for \code{\link{cv.rsv}}.  Setting
#' \code{n_v > 0} adds an overlap subsample (identifiable as
#' \code{S_e = 1 & S_o = 1}) that enables direct tests via
#' \code{\link{test_stability}} and \code{\link{test_no_direct_effect}}.
#'
#' @seealso \code{\link{cv.rsv}}, \code{\link{test_stability}},
#'   \code{\link{test_no_direct_effect}}
#'
#' @examples
#' # Main use case: Y unobserved in experimental sample
#' dat <- sim_rsv_data(n_e = 300, n_o = 700, seed = 42)
#' table(S_e = dat$S_e, S_o = dat$S_o, Y_observed = !is.na(dat$Y))
#'
#' # With overlap subsample for assumption testing
#' dat_v <- sim_rsv_data(n_e = 300, n_o = 700, n_v = 100, seed = 42)
#' table(S_e = dat_v$S_e, S_o = dat_v$S_o)
#'
#' @export
sim_rsv_data <- function(n_e = 300, n_o = 700, n_v = 0, tau = 0.10, p0 = 0.30,
                          n_r = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  stopifnot(
    is.numeric(n_e), n_e >= 1,
    is.numeric(n_o), n_o >= 1,
    is.numeric(n_v), n_v >= 0, n_v <= n_e,
    is.numeric(tau),
    is.numeric(p0), p0 > 0, p0 < 1,
    p0 + tau <= 1, p0 + tau >= 0
  )

  n_e     <- as.integer(n_e)
  n_o     <- as.integer(n_o)
  n_v     <- as.integer(n_v)
  n_pure  <- n_e - n_v          # pure experimental (S_e=1, S_o=0)
  n_total <- n_e + n_o

  # Row layout: [pure experimental | overlap | pure observational]
  S_e <- c(rep(1L, n_e), rep(0L, n_o))
  S_o <- c(rep(0L, n_pure), rep(1L, n_v), rep(1L, n_o))

  # Treatment: all experimental units (pure + overlap)
  D <- rep(NA_integer_, n_total)
  D[S_e == 1L] <- stats::rbinom(n_e, 1L, 0.5)

  # True Y for all units (latent for pure experimental)
  Y_true <- rep(NA_integer_, n_total)

  # Pure observational and overlap: Y observed
  obs_idx <- which(S_o == 1L & S_e == 0L)   # pure observational
  ovl_idx <- which(S_e == 1L & S_o == 1L)   # overlap
  pex_idx <- which(S_e == 1L & S_o == 0L)   # pure experimental

  Y_true[obs_idx] <- stats::rbinom(length(obs_idx), 1L, p0)
  Y_true[ovl_idx] <- stats::rbinom(length(ovl_idx), 1L,
                                    ifelse(D[ovl_idx] == 1L, p0 + tau, p0))
  Y_true[pex_idx] <- stats::rbinom(length(pex_idx), 1L,
                                    ifelse(D[pex_idx] == 1L, p0 + tau, p0))

  # RSVs generated from true Y (including latent experimental outcomes)
  R <- matrix(
    stats::rnorm(n_total * n_r, mean = 0.5 * Y_true, sd = 1),
    nrow = n_total, ncol = n_r
  )
  colnames(R) <- paste0("R", seq_len(n_r))

  # Observed Y: NA for pure experimental units
  Y_obs <- Y_true
  Y_obs[pex_idx] <- NA_integer_

  data.frame(Y = Y_obs, D = D, S_e = S_e, S_o = S_o, R)
}
