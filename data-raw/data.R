#' Andhra Pradesh Smartcard Study Data
#'
#' Base dataset from Muralidharan et al.'s smartcard study in Andhra Pradesh,
#' India, merged with SECC socioeconomic data. This dataset contains all 
#' villages with treatment assignment and outcomes, but excludes remote 
#' sensing features to reduce file size.
#'
#' @format A data frame with 8,312 rows and 10 columns:
#' \describe{
#'   \item{shrid2}{SHRUG village identifier (unique)}
#'   \item{spillover_20km}{Logical indicator for spillover-affected villages
#'     (within 20km of villages with different treatment status). Identified 
#'     using maximum bipartite matching and Konig's theorem to find the maximum 
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
#'   \item{Ylowinc}{Binary outcome: 1 if no households earn less than 5,000 
#'     rupees, 0 otherwise}
#'   \item{Ymidinc}{Binary outcome: 1 if no households earn more than 10,000 
#'     rupees, 0 otherwise}
#'   \item{clusters}{Cluster identifier (subdistrict name + district name in 
#'     Andhra Pradesh)}
#' }
#'
#' @details
#' This dataset must be merged with \code{remote_vars_p1} and 
#' \code{remote_vars_p2} to obtain the complete data with remote sensing 
#' features required for RSV estimation. The remote sensing variables are 
#' stored separately due to file size limits.
#'
#' @source
#' \itemize{
#'   \item Study data: Muralidharan, K., Niehaus, P., & Sukhtankar, S. (2016).
#'     Building State Capacity: Evidence from Biometric Smartcards in India.
#'     \emph{American Economic Review}, 106(10), 2895-2929. 
#'     \doi{10.1257/aer.20141346}
#'   \item SHRUG location data: Asher, S., Lunt, T., Matsuura, R., & Novosad, P. (2021).
#'     \url{https://dataverse.harvard.edu/api/access/datafile/10742739}
#'   \item SECC Consumption: 
#'     \url{https://dataverse.harvard.edu/api/access/datafile/10742743}
#'   \item SECC Income: 
#'     \url{https://dataverse.harvard.edu/api/access/datafile/10742876}
#' }
#'
#' @seealso 
#' \code{\link{remote_vars_p1}}, \code{\link{remote_vars_p2}} for remote 
#' sensing features
#'
#' @examples
#' # Load all three datasets
#' data("smartcard_data", package = "remoteoutcome")
#' data("remote_vars_p1", package = "remoteoutcome")
#' data("remote_vars_p2", package = "remoteoutcome")
#' force(smartcard_data)
#' force(remote_vars_p1)
#' force(remote_vars_p2)
#' 
#' # Merge to create complete dataset
#' smartcard_data <- smartcard_data %>%
#'   inner_join(remote_vars_p1, by = "shrid2") %>%
#'   inner_join(remote_vars_p2, by = "shrid2")
#'   
#' # Summary statistics
#' table(smartcard_data$D, useNA = "ifany")
#' table(smartcard_data$`Sample (Smartcard)`)
#' 
#' \dontrun{
#' # Use with RSV estimation
#' Y <- smartcard_data$Ycons
#' D <- smartcard_data$D
#' R <- smartcard_data %>% 
#'   select(starts_with("luminosity_"), starts_with("satellite_"))
#' S_e <- !is.na(D) & (rowSums(is.na(R)) == 0)
#' S_o <- !is.na(Y) & (rowSums(is.na(R)) == 0)
#' clusters <- smartcard_data$clusters
#' 
#' result <- rsv_estimate(
#'   Y = Y, D = D, S_e = S_e, S_o = S_o, R = R,
#'   method = "split",
#'   se_params = list(clusters = clusters)
#' )
#' }
"smartcard_data"


#' Remote Sensing Variables: VIIRS and MOSAIKS Features (Part 1)
#'
#' First part of remote sensing features for villages in the Andhra Pradesh 
#' smartcard study. Contains VIIRS nighttime lights data (2012-2021) and the 
#' first 2,000 MOSAIKS satellite imagery features.
#'
#' @format A data frame with 8,312 rows and 2,051 columns:
#' \describe{
#'   \item{shrid2}{SHRUG village identifier (unique) - merge key}
#'   \item{luminosity_min.*}{Minimum VIIRS nighttime luminosity (2012-2021)}
#'   \item{luminosity_max.*}{Maximum VIIRS nighttime luminosity (2012-2021)}
#'   \item{luminosity_mean.*}{Mean VIIRS nighttime luminosity (2012-2021)}
#'   \item{luminosity_sum.*}{Sum of VIIRS nighttime luminosity (2012-2021)}
#'   \item{luminosity_num_cells.*}{Number of VIIRS grid cells per village (2012-2021)}
#'   \item{satellite_*}{First 2,000 MOSAIKS satellite imagery features: 
#'     convolutional neural network features extracted from Planet imagery (2019)}
#' }
#'
#' @details
#' This dataset contains 50 VIIRS nighttime lights features (5 statistics × 
#' 10 years) and the first 2,000 of 4,000 MOSAIKS features. It must be merged 
#' with \code{smartcard_data} and \code{remote_vars_p2} to obtain the complete 
#' dataset for RSV estimation.
#'
#' The data is split into two parts (\code{remote_vars_p1} and 
#' \code{remote_vars_p2}) due to size limits while still 
#' providing access to all 4,000 MOSAIKS features.
#'
#' @source
#' \itemize{
#'   \item VIIRS nighttime lights: 
#'     \url{https://dataverse.harvard.edu/api/access/datafile/10742856}
#'   \item MOSAIKS features: Rolf, E., Proctor, J., Carleton, T., Bolliger, I., 
#'     Shankar, V., Ishihara, M., Recht, B., & Hsiang, S. (2021). 
#'     A generalizable and accessible approach to machine learning with global 
#'     satellite imagery. \emph{Nature Communications}, 12, 4392.
#'     \doi{10.1038/s41467-021-24638-z}
#' }
#'
#' @seealso 
#' \code{\link{smartcard_data}} for base data,
#' \code{\link{remote_vars_p2}} for remaining MOSAIKS features,
#' \code{vignette("construct-remote-vars", package = "remoteoutcome")} for 
#' instructions on generating these features
#'
#' @examples
#' # Load all datasets
#' data("smartcard_data", package = "remoteoutcome")
#' data("remote_vars_p1", package = "remoteoutcome")
#' data("remote_vars_p2", package = "remoteoutcome")
#' force(smartcard_data)
#' force(remote_vars_p1)
#' force(remote_vars_p2)
#' 
#' # Merge to create complete dataset
#' smartcard_data <- smartcard_data %>%
#'   inner_join(remote_vars_p1, by = "shrid2") %>%
#'   inner_join(remote_vars_p2, by = "shrid2")
#'   
#' # Check dimensions
#' dim(smartcard_data)  # Should be 8,312 × 4,060
#' 
#' # Summary of VIIRS features
#' summary(remote_vars_p1[, grep("luminosity_mean", names(remote_vars_p1))])
#' 
#' # Summary of first 10 satellite features
#' summary(remote_vars_p1[, paste0("satellite_", 1:10)])
"remote_vars_p1"


#' Remote Sensing Variables: MOSAIKS Features (Part 2)
#'
#' Second part of remote sensing features for villages in the Andhra Pradesh 
#' smartcard study. Contains the remaining 2,000 MOSAIKS satellite imagery 
#' features (features 2001-4000).
#'
#' @format A data frame with 8,312 rows and 2,001 columns:
#' \describe{
#'   \item{shrid2}{SHRUG village identifier (unique) - merge key}
#'   \item{satellite_2001 to satellite_4000}{Remaining 2,000 MOSAIKS satellite 
#'     imagery features: convolutional neural network features extracted from 
#'     Planet imagery (2019)}
#' }
#'
#' @details
#' This dataset contains the second half of the 4,000 MOSAIKS features. It 
#' must be merged with \code{smartcard_data} and \code{remote_vars_p1} to 
#' obtain the complete dataset for RSV estimation.
#'
#' The MOSAIKS features are high-dimensional representations of satellite 
#' imagery that capture various aspects of land use, development, vegetation, 
#' and built environment. These features are derived from applying 
#' convolutional neural networks to Planet imagery from 2019.
#'
#' @source
#' MOSAIKS features: Rolf, E., Proctor, J., Carleton, T., Bolliger, I., 
#' Shankar, V., Ishihara, M., Recht, B., & Hsiang, S. (2021). 
#' A generalizable and accessible approach to machine learning with global 
#' satellite imagery. \emph{Nature Communications}, 12, 4392.
#' \doi{10.1038/s41467-021-24638-z}
#'
#' @seealso 
#' \code{\link{smartcard_data}} for base data,
#' \code{\link{remote_vars_p1}} for VIIRS features and first 2000 MOSAIKS features,
#' \code{vignette("construct-remote-vars", package = "remoteoutcome")} for 
#' instructions on generating these features
#'
#' @examples
#' # Load all datasets
#' data("smartcard_data", package = "remoteoutcome")
#' data("remote_vars_p1", package = "remoteoutcome")
#' data("remote_vars_p2", package = "remoteoutcome")
#' force(smartcard_data)
#' force(remote_vars_p1)
#' force(remote_vars_p2)
#' 
#' # Merge to create complete dataset
#' smartcard_data <- smartcard_data %>%
#'   inner_join(remote_vars_p1, by = "shrid2") %>%
#'   inner_join(remote_vars_p2, by = "shrid2")
#'   
#' # Verify all 4,000 satellite features are present
#' satellite_cols <- grep("^satellite_", names(smartcard_data), value = TRUE)
#' length(satellite_cols)  # Should be 4000
#' 
#' # Summary of last 10 satellite features
#' summary(remote_vars_p2[, paste0("satellite_", 3991:4000)])
"remote_vars_p2"



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
#'   \item \code{\link{smartcard_data}} for the underlying data without predictions
#'   \item \code{vignette(package = "remoteoutcome")} for examples
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