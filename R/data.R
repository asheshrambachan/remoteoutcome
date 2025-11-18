#' Predicted Values and Observations for Consumption Outcome (Real Sample Split)
#'
#' Pre-computed predictions and observations from the RSV estimation procedure
#' for the consumption outcome (Ycons) using the real experimental/observational
#' sample split. This dataset allows users to run \code{\link{rsv_estimate}} with
#' pre-fitted predictions, avoiding the computational cost of re-fitting random
#' forest models.
#'
#' @format A data frame with 8,312 rows and 9 columns:
#' \describe{
#'   \item{Y}{Binary outcome: consumption indicator (1 = bottom quartile, 0 = otherwise).
#'     \code{NA} for experimental-only observations where outcomes are not observed.}
#'   \item{D}{Binary treatment indicator (1 = treated, 0 = control).
#'     \code{NA} for observational-only observations where treatment is not assigned.}
#'   \item{S_e}{Experimental sample indicator: 1 if unit is in experimental sample
#'     (treatment observed), 0 otherwise.}
#'   \item{S_o}{Observational sample indicator: 1 if unit is in observational sample
#'     (outcome observed), 0 otherwise.}
#'   \item{pred_Y}{Predicted probability \eqn{P(Y = 1 \mid R, S_o = 1)}, where
#'     \eqn{R} includes VIIRS nighttime lights and MOSAIKS satellite features.
#'     Fitted using random forest on observational sample.}
#'   \item{pred_D}{Predicted probability \eqn{P(D = 1 \mid R, S_e = 1)}.
#'     Fitted using random forest on experimental sample.}
#'   \item{pred_S_e}{Predicted probability \eqn{P(S_e = 1 \mid R)}.
#'     Fitted using random forest on full sample.}
#'   \item{pred_S_o}{Predicted probability \eqn{P(S_o = 1 \mid R)}.
#'     Fitted using random forest on full sample.}
#'   \item{clusters}{Cluster identifier (subdistrict name + district name) for
#'     computing cluster-robust standard errors.}
#' }
#'
#' @details
#' \strong{Generation Process:}
#' 
#' This dataset was generated using the following procedure:
#' \enumerate{
#'   \item \strong{Data preparation:}
#'     \itemize{
#'       \item Started with \code{data_real} (real experimental/observational split)
#'       \item Merged with satellite features (MOSAIKS 4,000-dimensional features)
#'       \item Constructed remotely sensed variable \eqn{R} from VIIRS nighttime
#'         lights (2012-2021) and MOSAIKS features
#'     }
#'   \item \strong{Sample indicators:}
#'     \itemize{
#'       \item \code{S_e = 1} for units where treatment \code{D} and covariates
#'         \code{R} are both observed (experimental sample)
#'       \item \code{S_o = 1} for units where outcome \code{Y} and covariates
#'         \code{R} are both observed (observational sample)
#'       \item Some units have both \code{S_e = 1} and \code{S_o = 1} (overlap sample)
#'     }
#'   \item \strong{Prediction fitting:}
#'     \itemize{
#'       \item Used \code{method = "none"} (no sample splitting - all data used
#'         for both training and prediction)
#'       \item Fitted four random forest models using \code{ranger}:
#'         \itemize{
#'           \item \code{pred_Y}: \eqn{P(Y = 1 \mid R, S_o = 1)} fitted on observational sample
#'           \item \code{pred_D}: \eqn{P(D = 1 \mid R, S_e = 1)} fitted on experimental sample
#'           \item \code{pred_S_e}: \eqn{P(S_e = 1 \mid R)} fitted on full sample
#'           \item \code{pred_S_o}: \eqn{P(S_o = 1 \mid R)} fitted on full sample
#'         }
#'       \item Random forest parameters: 100 trees, class weights \code{c(10, 1)}
#'         for \code{pred_Y} model (upweighting rare outcome), seed = 42
#'     }
#'   \item \strong{Initial estimate:}
#'     \itemize{
#'       \item Computed \code{theta_init} (initial treatment effect estimate) on
#'         the training data using the predictions
#'       \item Stored as an attribute: \code{attr(pred_real_Ycons, "theta_init")}
#'     }
#' }
#'
#'
#' @section Attribute:
#' The dataset has one attribute:
#' \describe{
#'   \item{theta_init}{Initial estimate of the treatment effect computed on the
#'     training data. This is used as a starting value in the RSV estimation
#'     procedure. Access with \code{attr(pred_real_Ycons, "theta_init")}.}
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{rsv_estimate}} for using this dataset to estimate treatment effects
#' }
#'
#' @examples
#' # Load the dataset
#' data("pred_real_Ycons", package = "remoteoutcome")
#'
#' # Examine structure
#' str(pred_real_Ycons)
#'
#' # Check sample sizes
#' table(S_e = pred_real_Ycons$S_e, S_o = pred_real_Ycons$S_o)
#'
#' # View prediction distributions
#' summary(pred_real_Ycons[, c("pred_Y", "pred_D", "pred_S_e", "pred_S_o")])
#'
#' # Access the initial treatment effect estimate
#' attr(pred_real_Ycons, "theta_init")
#'
#' \dontrun{
#' # Estimate treatment effect using pre-computed predictions
#' result <- rsv_estimate(
#'   Y = pred_real_Ycons$Y,
#'   D = pred_real_Ycons$D,
#'   S_e = pred_real_Ycons$S_e,
#'   S_o = pred_real_Ycons$S_o,
#'   pred_Y = pred_real_Ycons$pred_Y,
#'   pred_D = pred_real_Ycons$pred_D,
#'   pred_S_e = pred_real_Ycons$pred_S_e,
#'   pred_S_o = pred_real_Ycons$pred_S_o,
#'   method = "predictions",
#'   theta_init = attr(pred_real_Ycons, "theta_init"),
#'   se = TRUE,
#'   se_params = list(
#'     B = 1000,
#'     fix_seed = TRUE,
#'     clusters = pred_real_Ycons$clusters
#'   ),
#'   cores = 7
#' )
#'
#' # View results
#' print(result)
#' confint(result)
#' }
#'
#' @keywords datasets
"pred_real_Ycons"