library(dplyr)
library(remoteoutcome)

# Load the data
data("smartcard_data_p1", package="remoteoutcome")
data("smartcard_data_p2", package="remoteoutcome")

# Merge to create complete dataset
smartcard_data <- inner_join(smartcard_data_p1, smartcard_data_p2, by="shrid2")
rm(smartcard_data_p1, smartcard_data_p2)

data_real <- create_data_real(smartcard_data)

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
usethis::use_data(pred_real_Ycons, overwrite = TRUE, compress = "xz")
