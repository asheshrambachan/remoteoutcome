#' Andhra Pradesh Smartcard Study Data
#'
#' Complete dataset from Muralidharan et al.'s smartcard study in Andhra Pradesh,
#' India, merged with SECC socioeconomic data. 
#'
#' @format A data frame with 8,312 rows and 10 columns:
#' \describe{
#'   \item{shrid2}{SHRUG village identifier (unique)}
#'   \item{spillover_20km}{Logical indicator for spillover-affected villages
#'     (within 20km of villages with different treatment status). The identification 
#'     uses maximum bipartite matching and Konig's theorem to find the maximum 
#'     independent set.}
#'   \item{tot_p}{Total population from SECC census}
#'   \item{tot_f}{Total number of families from SECC census}
#'   \item{Sample (Smartcard)}{Factor indicating sample assignment:
#'     "Experimental: Treated (2010)", "Experimental: Untreated (2011)",
#'     "Experimental: Untreated (2012)", or "Observational (N/A)"}
#'   \item{D}{Treatment indicator: 
#'     1 if "Experimental: Treated (2010)" wave, 
#'     0 if "Experimental: Untreated (2011)" or "Experimental: Untreated (2012)" waves, 
#'     NA if "Observational (N/A)" wave}
#'   \item{Ycons}{Binary outcome: 1 if village consumption is in bottom quartile
#'     (<= 18,946.61 rupees per capita), 0 otherwise}
#'   \item{Ylowinc}{Binary outcome: 1 if no households earn less than 5,000 rupees, 0 otherwise}
#'   \item{Ymidinc}{Binary outcome: 1 if no households earn more than 10,000 rupees, 0 otherwise}
#'   \item{clusters}{Cluster identifier (Subdistrict (mandal) name + District name in Andhra Pradesh)}
#' }
#'
#' @source
#' \itemize{
#'   \item Study data: Muralidharan, K., Niehaus, P., & Sukhtankar, S. (2016).
#'     Building State Capacity: Evidence from Biometric Smartcards in India.
#'     \emph{American Economic Review}, 106(10), 2895-2929. 
#'     \url{https://www.openicpsr.org/openicpsr/project/113012/version/V1/view?path=/openicpsr/113012/fcr:versions/V1/20141346_data/data/balance-for-ap-mandal-comparison.dta&type=file}
#'   \item SHRUG location data: \url{https://dataverse.harvard.edu/api/access/datafile/10742739}
#'   \item Socio-Economic and Caste Census (SECC) Consumption Data: \url{https://dataverse.harvard.edu/api/access/datafile/10742743}
#'   \item Socio-Economic and Caste Census (SECC) Income Data: \url{https://dataverse.harvard.edu/api/access/datafile/10742876}
#'   \item VIIRS: \url{https://dataverse.harvard.edu/api/access/datafile/10742856}
#'   \item SHRUG shapefiles: Asher, S., Lunt, T., Matsuura, R., & Novosad, P. (2021).
#'     Development research at high geographic resolution: An analysis of night-lights,
#'     firms, and poverty in India using the SHRUG open data platform.
#'     \emph{The World Bank Economic Review}, 35(4), 845-871.
#'     \url{https://www.devdatalab.org/shrug_download}
#'   \item MOSAIKS: Rolf, E., Proctor, J., Carleton, T., Bolliger, I., Shankar, V.,
#'     Ishihara, M., Recht, B., & Hsiang, S. (2021). A generalizable and accessible
#'     approach to machine learning with global satellite imagery.
#'     \emph{Nature Communications}, 12, 4392.
#'     \doi{10.1038/s41467-021-24638-z}
#' }
#'
#' @examples
#' \dontrun{
#' # Load the data
#' data(smartcard_data, library="remoteoutcome")
#' data(remote_vars_p1, library="remoteoutcome")
#' data(remote_vars_p2, library="remoteoutcome")
#' 
#' smartcard_data <- smartcard_data %>%
#'   inner_join(remote_vars_p1, by="shrid2") %>%
#'   inner_join(remote_vars_p2, by="shrid2")
#'   
#' # Summary of treatment assignment
#' table(smartcard_data$D, useNA = "ifany")
#'
#' # Villages affected by spillovers
#' table(smartcard_data$spillover_20km)
#'
#' # Sample distribution
#' table(smartcard_data$`Sample (Smartcard)`)
#' }
"smartcard_data"


#' Andhra Pradesh Smartcard Study Data with Remote Sensing Features
#'
#' Complete dataset from Muralidharan et al.'s smartcard study in Andhra Pradesh,
#' India, merged with SECC socioeconomic data and VIIRS nighttime lights data. 
#' This dataset contains all villages with treatment assignment, and outcomes.
#'
#' @format A data frame with 8,312 rows and 60 columns:
#' \describe{
#'   \item{shrid2}{SHRUG village identifier (unique)}
#'   \item{satellite_*}{MOSAIKS satellite imagery features: 4,000-dimensional  
#'     convolutional neural network features extracted from Planet imagery (2019)}
#' }
#'
#' @source
#' \itemize{
#'   \item Study data: Muralidharan, K., Niehaus, P., & Sukhtankar, S. (2016).
#'     Building State Capacity: Evidence from Biometric Smartcards in India.
#'     \emph{American Economic Review}, 106(10), 2895-2929. 
#'     \url{https://www.openicpsr.org/openicpsr/project/113012/version/V1/view?path=/openicpsr/113012/fcr:versions/V1/20141346_data/data/balance-for-ap-mandal-comparison.dta&type=file}
#'   \item SHRUG location data: \url{https://dataverse.harvard.edu/api/access/datafile/10742739}
#'   \item Socio-Economic and Caste Census (SECC) Consumption Data: \url{https://dataverse.harvard.edu/api/access/datafile/10742743}
#'   \item Socio-Economic and Caste Census (SECC) Income Data: \url{https://dataverse.harvard.edu/api/access/datafile/10742876}
#'   \item VIIRS: \url{https://dataverse.harvard.edu/api/access/datafile/10742856}
#'   \item SHRUG shapefiles: Asher, S., Lunt, T., Matsuura, R., & Novosad, P. (2021).
#'     Development research at high geographic resolution: An analysis of night-lights,
#'     firms, and poverty in India using the SHRUG open data platform.
#'     \emph{The World Bank Economic Review}, 35(4), 845-871.
#'     \url{https://www.devdatalab.org/shrug_download}
#'   \item MOSAIKS: Rolf, E., Proctor, J., Carleton, T., Bolliger, I., Shankar, V.,
#'     Ishihara, M., Recht, B., & Hsiang, S. (2021). A generalizable and accessible
#'     approach to machine learning with global satellite imagery.
#'     \emph{Nature Communications}, 12, 4392.
#'     \doi{10.1038/s41467-021-24638-z}
#' }
#'
#' @examples
#' \dontrun{
#' # Load the data
#' data(smartcard_data, library="remoteoutcome")
#' data(remote_vars_p1, library="remoteoutcome")
#' data(remote_vars_p2, library="remoteoutcome")
#' data(remote_vars_p3, library="remoteoutcome")
#' 
#' smartcard_data <- smartcard_data %>%
#'   inner_join(remote_vars_p1, by="shrid2") %>%
#'   inner_join(remote_vars_p2, by="shrid2") %>%
#'   inner_join(remote_vars_p3, by="shrid2")
#'   
#' # Summary of treatment assignment
#' table(smartcard_data$D, useNA = "ifany")
#'
#' # Villages affected by spillovers
#' table(smartcard_data$spillover_20km)
#'
#' # Sample distribution
#' table(smartcard_data$`Sample (Smartcard)`)
#' }
"remote_vars_p3"


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
#'           \item \code{pred_Y}: \eqn{E[Y \mid R, S_o = 1]} fitted on observational sample
#'           \item \code{pred_D}: \eqn{E[D \mid R, S_e = 1]} fitted on experimental sample
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
#'   \item \code{\link{smartcard_data}} for the underlying data without predictions
#'   \item \code{vignette(package = "remoteoutcome")} for examples
#' }
#'
#' @examples
#' # Load the dataset
#' data(pred_real_Ycons, package = "remoteoutcome")
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