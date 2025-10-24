# =============================================================================
# Common Utilities
# =============================================================================
#' 
#' #' Generate synthetic RSV data
#' #'
#' #' Generates synthetic data for testing and simulating RSV estimation methods.
#' #' The data generating process creates treatment assignment, outcomes, and
#' #' remotely sensed variables under the RSV framework.
#' #'
#' #' @param n Sample size.
#' #' @param tau Treatment effect parameter.
#' #' @param X Covariate matrix.
#' #' @param D Treatment indicator vector.
#' #' @param Y Outcome vector.
#' #'
#' #' @return A list containing:
#' #' \describe{
#' #'   \item{X}{Covariate matrix}
#' #'   \item{D}{Treatment indicator (0/1)}
#' #'   \item{Y}{Outcome variable (NA for experimental sample)}
#' #'   \item{Se}{Logical indicator for experimental sample}
#' #'   \item{So}{Logical indicator for observational sample}
#' #' }
#' #'
#' #' @export
#' #'
#' #' @examples
#' #' # Generate synthetic data
#' #' n <- 1000
#' #' X <- matrix(rnorm(n * 5), ncol = 5)
#' #' D <- rbinom(n, 1, 0.5)
#' #' Y <- rbinom(n, 1, 0.3)
#' #' data <- generate_rsv_data(n, tau = 0.1, X, D, Y)
#' generate_rsv_data <- function(n, tau, X, D, Y) {
#'   Y1_id <- which(Y == 1)
#'   Y0_id <- which(Y == 0)
#'   p_Y1_D1 <- mean(Y[D == 1], na.rm = TRUE) + tau
#'   p_Y1_D0 <- mean(Y[D == 0], na.rm = TRUE)
#' 
#'   # Draw new D and outcomes
#'   D <- rbinom(n, 1, prob = 0.5)
#'   Y_D1 <- rbinom(n, 1, prob = p_Y1_D1)
#'   Y_D0 <- rbinom(n, 1, prob = p_Y1_D0)
#'   X_Y1_id <- sample(Y1_id, size = n, replace = TRUE)
#'   X_Y0_id <- sample(Y0_id, size = n, replace = TRUE)
#' 
#'   # Define Y
#'   Y <- D * Y_D1 + (1 - D) * Y_D0
#' 
#'   # Define X
#'   X_Y1 <- X[X_Y1_id, , drop = FALSE]
#'   X_Y0 <- X[X_Y0_id, , drop = FALSE]
#'   X <- Y * X_Y1 + (1 - Y) * X_Y0
#' 
#'   # Define Se and So
#'   Se <- rep(TRUE, n)  # S = e for D = 0 and D = 1
#'   So <- D != 1        # S = o for D = 0
#'   Y <- ifelse(So == TRUE, Y, NA)
#' 
#'   list(X = X, D = D, Y = Y, Se = Se, So = So)
#' }

# #' Create experimental/observational indicators and mask variables
# #'
# #' Converts a dataset with a smartcard sample indicator into:
# #' - `S`: sample membership ("e", "o", "both")
# #' - `D`: treatment observed only for experimental/both samples
# #' - outcomes observed only for observational/both samples
# #'
# #' @param data A data.frame  from 
# #'
# #' @return The input `data` with added/modified columns `S`, `D`, and outcomes.
# #' 
# #' @export
# #' @importFrom dplyr mutate case_when
# #' @importFrom rlang .data
# generate_real_data <- function(data){
#   data %>%
#   dplyr::mutate(
#     # Sample indicator
#     S = dplyr::case_when(
#       `Sample (Smartcard)` == "Experimental: Treated (2010)" ~ "e", # Experimental only
#       `Sample (Smartcard)` == "Experimental: Untreated (2011)" ~ "both",  # Both samples
#       `Sample (Smartcard)` == "Experimental: Untreated (2012)" ~ "e", # Experimental only
#       `Sample (Smartcard)` == "Observational (N/A)" ~ "o" # Observational only
#     ),
    
#     # Treatment observed only in experimental sample
#     D = ifelse((S=="both") | (S=="e"), D, NA),
    
#     # Outcomes observed only in observational sample
#     Ycons = ifelse((S=="both") | (S=="o"), Ycons, NA),
#     Ylowinc = ifelse((S=="both") | (S=="o"), Ylowinc, NA),
#     Ymidinc = ifelse((S=="both") | (S=="o"), Ymidinc, NA)
#   )
# }

# #' Create experimental/observational indicators and mask variables
# #'
# #' Converts a dataset with a smartcard sample indicator into:
# #' - `S`: sample membership ("e", "o", "both")
# #' - `D`: treatment observed only for experimental/both samples
# #' - outcomes observed only for observational/both samples
# #'
# #' @param data A data.frame  from 
# #'
# #' @return The input `data` with added/modified columns `S`, `D`, and outcomes.
# #' 
# #' @export
# #' @importFrom dplyr mutate case_when
# #' @importFrom rlang .data
# generate_synth_data <- function(data){
#   dplyr::mutate(
#     data,
    
#     # Sample indicators
#     S_e = `Sample (Smartcard)` %in% c(
#       "Experimental: Treated (2010)", 
#       "Experimental: Untreated (2011)", 
#       "Experimental: Untreated (2012)"
#     ), # Experimental sample (has D)
    
#     S_o = clusters %in% obs_clusters, # Synthetic observational sample

#     # Sample indicator for rsv_estimate (categorical)
#     S = case_when(
#       S_e & S_o ~ "both",    # Both samples
#       S_e & !S_o ~ "e",      # Experimental only
#       !S_e & S_o ~ "o"       # Observational only
#     ),
    
#     # Treatment observed only in experimental sample
#     D = ifelse((S=="both") | (S=="e"), D, NA),
    
#     # Outcomes observed only in observational sample
#     Ycons = ifelse((S=="both") | (S=="o"), Ycons, NA),
#     Ylowinc = ifelse((S=="both") | (S=="o"), Ylowinc, NA),
#     Ymidinc = ifelse((S=="both") | (S=="o"), Ymidinc, NA)
#   )
# }