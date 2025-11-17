# =============================================================================
# Data Constructor Functions
# Functions to create experimental/observational sample splits from smartcard_data
# =============================================================================

#' @importFrom dplyr %>% mutate case_when if_else select starts_with
NULL


#' Create Real Experimental/Observational Split Dataset
#'
#' Transforms \code{smartcard_data} into the real experimental/observational split 
#' where treatment is observed in "Experimental: Treated (2010)", "Experimental: 
#' Untreated (2011)", and "Experimental: Untreated (2012)" waves and outcomes are 
#' observed in the "Experimental: Untreated (2011)" and "Observational (N/A)" waves.
#'
#' @param data \code{smartcard_data} data frame included in the package.
#'
#' @return A data frame similar to \code{smartcard_data}, containing:
#' \describe{
#'   \item{shrid2}{SHRUG village identifier}
#'   \item{spillover_20km}{Spillover indicator}
#'   \item{S}{Sample indicator: "e" = experimental only, "o" = observational only,
#'     "both" = experimental and observational samples}
#'   \item{D}{Treatment indicator, \code{NA} when not observed (only observed 
#'     when S = "e" or "both")}
#'   \item{Ycons, Ylowinc, Ymidinc}{Outcome variables, \code{NA} when not 
#'     observed (only observed when S = "o" or "both")}
#'   \item{clusters}{Cluster identifier (Subdistrict (mandal) name + District name in Andhra Pradesh)}
#'   \item{luminosity_*}{VIIRS nighttime lights features}
#'   \item{satellite_*}{MOSAIKS satellite feature}
#' }
#' 
#' @details
#' The function implements the real experimental design through the following steps:
#' \enumerate{
#'   \item Villages from "Experimental: Treated (2010)", "Experimental: Untreated (2011)", 
#'     "Experimental: Untreated (2012)" waves are marked as S_e = TRUE
#'   \item Villages from "Experimental: Untreated (2011)" and "Observational (N/A)" 
#'     waves are marked as S_o = TRUE
#'   \item Final sample indicator is created:
#'     \itemize{
#'       \item S = "both" if S_e = TRUE and S_o = TRUE (experimental and observational villages)
#'       \item S = "e" if S_e = TRUE and S_o = FALSE (experimental only)
#'       \item S = "o" if S_e = FALSE and S_o = TRUE (observational only)
#'     }
#'   \item Treatment visible only when S = "e" or "both"; 
#'     outcomes visible only when S = "o" or "both"
#' }
#' 
#' @seealso \code{\link{smartcard_data}} for the input dataset structure,
#'   \code{\link{create_data_synth}} for creating synthetic sample splits
#'
#' @examples
#' \dontrun{
#' # Load complete dataset
#' data(smartcard_data, package = "remoteoutcome")
#'
#' # Create real experimental/observational split
#' data_real <- create_data_real(smartcard_data)
#'
#' # Sample membership distribution
#' table(data_real$S)
#'
#' # Treatment and outcome overlap
#' with(data_real, table(observe_D = !is.na(D), observe_Y = !is.na(Ycons)))
#' }
#'
#' @export
create_data_real <- function(data) {
 
  # Transform data
  result <- data %>%
    dplyr::mutate(
      # Sample indicator based on original sample assignment
      S = dplyr::case_when(
        `Sample (Smartcard)` == "Experimental: Treated (2010)" ~ "e",
        `Sample (Smartcard)` == "Experimental: Untreated (2011)" ~ "both",
        `Sample (Smartcard)` == "Experimental: Untreated (2012)" ~ "e",
        `Sample (Smartcard)` == "Observational (N/A)" ~ "o"
      ),
      
      # Treatment observed only when S = "e" or S = "both"
      D = dplyr::if_else(S %in% c("both", "e"), D, NA),
      
      # Outcomes observed only when S = "o" or S = "both"
      Ycons = dplyr::if_else(S %in% c("both", "o"), Ycons, NA),
      Ylowinc = dplyr::if_else(S %in% c("both", "o"), Ylowinc, NA),
      Ymidinc = dplyr::if_else(S %in% c("both", "o"), Ymidinc, NA)
    ) %>%
    dplyr::select(
      shrid2, spillover_20km, S, D, Ycons, Ylowinc, Ymidinc, clusters,
      dplyr::starts_with("luminosity_"), 
      dplyr::starts_with("satellite_")
    )
  
  # Validation checks
  n_exp <- sum(result$S == "e", na.rm = TRUE)
  n_obs <- sum(result$S == "o", na.rm = TRUE)
  n_both <- sum(result$S == "both", na.rm = TRUE)
  
  if (n_exp == 0 || n_obs == 0) {
    warning("Sample split resulted in empty experimental or observational sample")
  }
  
  cat(sprintf(
    "Created data_real with %d observations:\n  Experimental only: %d\n  Observational only: %d\n  Both samples: %d\n",
    nrow(result), n_exp, n_obs, n_both
  ))
  
  return(result)
}


#' Create Synthetic Experimental/Observational Split Dataset
#'
#' Transforms the complete dataset (\code{smartcard_data}) into a synthetic 
#' experimental/observational split where top 50% of the clusters are assigned to 
#' observational sample. This creates an artificial separation for testing 
#' RSV methods.
#'
#' @param data \code{smartcard_data} data frame included in the package.
#'   
#' @return A data frame similar to \code{smartcard_data}, containing:
#' \describe{
#'   \item{shrid2}{SHRUG village identifier}
#'   \item{spillover_20km}{Spillover indicator}
#'   \item{S}{Sample indicator: "e" = experimental only, "o" = observational only,
#'     "both" = experimental and observational samples}
#'   \item{D}{Treatment indicator, \code{NA} when not observed (only observed 
#'     when S = "e" or "both")}
#'   \item{Ycons, Ylowinc, Ymidinc}{Outcome variables, \code{NA} when not 
#'     observed (only observed when S = "o" or "both")}
#'   \item{clusters}{Cluster identifier (Subdistrict (mandal) name + District name in Andhra Pradesh)}
#'   \item{luminosity_*}{VIIRS nighttime lights features}
#'   \item{satellite_*}{MOSAIKS satellite features}
#' }
#'
#' @details
#' The function creates a synthetic experimental/observational split through 
#' the following steps:
#' \enumerate{
#'   \item Original experimental villages are marked as S_e = TRUE
#'   \item Top 50% of clusters are marked as S_o = TRUE
#'   \item Final sample indicator is created:
#'     \itemize{
#'       \item S = "both" if S_e = TRUE and S_o = TRUE (experimental and observational villages)
#'       \item S = "e" if S_e = TRUE and S_o = FALSE (experimental only)
#'       \item S = "o" if S_e = FALSE and S_o = TRUE (observational only)
#'     }
#'   \item Treatment visible only when S = "e" or "both"; 
#'     outcomes visible only when S = "o" or "both"
#' }
#'
#' @seealso \code{\link{smartcard_data}} for the input dataset structure,
#'   \code{\link{create_data_real}} for creating real experimental splits
#'
#' @examples
#' \dontrun{
#' # Load complete dataset
#' data(smartcard_data, package = "remoteoutcome")
#'
#' data_synth <- create_data_synth(smartcard_data)
#'
#' # Sample distribution
#' table(data_synth$S)
#' 
#' # Treatment and outcome overlap
#' with(data_synth, table(observe_D = !is.na(D), observe_Y = !is.na(Ycons)))
#' }
#'
#' @export
create_data_synth <- function(data) {

  # Assign half the clusters to observational sample
  clusters <- unique(data$clusters)
  obs_clusters <- clusters[1:floor(length(clusters) / 2)]
  n_obs_clusters <- length(obs_clusters)
  
  # Transform data
  result <- data %>%
    dplyr::mutate(
      # Experimental sample: original experimental waves (1, 2, 3)
      S_e = `Sample (Smartcard)` %in% c(
        "Experimental: Treated (2010)", 
        "Experimental: Untreated (2011)", 
        "Experimental: Untreated (2012)"
      ),
      
      # Observational sample: randomly selected clusters
      S_o = clusters %in% obs_clusters,
      
      # Combined sample indicator
      S = dplyr::case_when(
        S_e & S_o ~ "both",
        S_e & !S_o ~ "e",
        !S_e & S_o ~ "o"
      ),
      
      # Treatment observed only when S = "e" or S = "both"
      D = dplyr::if_else(S %in% c("both", "e"), D, NA),
      
      # Outcomes observed only when S = "o" or S = "both"
      Ycons = dplyr::if_else(S %in% c("both", "o"), Ycons, NA),
      Ylowinc = dplyr::if_else(S %in% c("both", "o"), Ylowinc, NA),
      Ymidinc = dplyr::if_else(S %in% c("both", "o"), Ymidinc, NA)
    ) %>%
    dplyr::select(
      shrid2, spillover_20km, S, D, Ycons, Ylowinc, Ymidinc, clusters,
      dplyr::starts_with("luminosity_"), 
      dplyr::starts_with("satellite_")
    )
  
  # Validation checks
  n_exp <- sum(result$S == "e", na.rm = TRUE)
  n_obs <- sum(result$S == "o", na.rm = TRUE)
  n_both <- sum(result$S == "both", na.rm = TRUE)
  
  if (n_exp == 0 || n_obs == 0) {
    warning("Sample split resulted in empty experimental or observational sample")
  }
  
  cat(sprintf(
    "Created data_synth with %d observations:\n  Experimental only: %d\n  Observational only: %d\n  Both samples: %d\n  Clusters in observational sample: %d of %d (%.1f%%)\n",
    nrow(result), n_exp, n_obs, n_both, n_obs_clusters, length(clusters), 
    n_obs_clusters/length(clusters) * 100
  ))
  
  return(result)
}
