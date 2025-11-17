# ==============================================================================
# USER CONFIGURATION - EDIT THIS PATHS
# ==============================================================================

path_to_satellite_features <- "~/Documents/remoteoutcome/data-raw/satellite_features.rds"

# ==============================================================================
# END USER CONFIGURATION
# ==============================================================================

library(dplyr)
library(remoteoutcome)

satellite_features <- readRDS(path_to_satellite_features) %>%
  select(-contains("lat"), -contains("lon"))

data("data_real", package = "remoteoutcome")
data_real <- data_real %>%
  inner_join(satellite_features, by = "shrid2")

Y <- data_real$Ycons
D <- data_real$D
R <- data_real %>% select(starts_with("luminosity"), starts_with("satellite"))
S_e <- !is.na(D) & (rowSums(is.na(R)) == 0)
S_o <- !is.na(Y) & (rowSums(is.na(R)) == 0)
clusters <- data_real$clusters

rsv_real_Ycons <- rsv_estimate(
  Y = Y,
  D = D,
  R = R,
  S_e = !is.na(D) & (rowSums(is.na(R)) == 0),
  S_o = !is.na(Y) & (rowSums(is.na(R)) == 0),
  eps = 1e-2,
  method = "none",
  ml_params = list(ntree = 100, classwt_Y = c(10, 1), seed=42),
  se = TRUE,
  se_params = list(B = 1000, fix_seed = TRUE, clusters = clusters),
  cores = cores
)

pred_real_Ycons <- data.frame(
  Y = rsv_real_Ycons$observations$Y,
  D = rsv_real_Ycons$observations$D,
  S_e = rsv_real_Ycons$observations$S_e,
  S_o = rsv_real_Ycons$observations$S_o,
  pred_Y = rsv_real_Ycons$predictions$Y,
  pred_D = rsv_real_Ycons$predictions$D,
  pred_S_e = rsv_real_Ycons$predictions$S_e,
  pred_S_o = rsv_real_Ycons$predictions$S_o,
  clusters = rsv_real_Ycons$clusters
)
attr(pred_real_Ycons, "theta_init") <- rsv_real_Ycons$theta_init
# usethis::use_data(pred_real_Ycons, overwrite = TRUE)
save(pred_real_Ycons, file = "data/pred_real_Ycons.rda", compress = "xz", compression_level = 9)