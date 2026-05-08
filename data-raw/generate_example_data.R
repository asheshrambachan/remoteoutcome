## Script to regenerate the rsv_example dataset bundled with the package.
## Run this script from the package root with:
##   source("data-raw/generate_example_data.R")

source("R/sim_rsv_data.R")   # standalone, no other package R/ dependencies

rsv_example <- sim_rsv_data(
  n_e  = 300,
  n_o  = 700,
  n_v  = 100,
  tau  = 0.10,
  p0   = 0.30,
  n_r  = 5,
  seed = 42
)

usethis::use_data(rsv_example, overwrite = TRUE)
