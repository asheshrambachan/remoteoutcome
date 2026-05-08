# Overlap design: 150 experimental units also in observational sample (S_e=1, S_o=1)
dat <- sim_rsv_data(n_e = 300, n_o = 700, n_v = 150, seed = 42)
R   <- as.matrix(dat[, paste0("R", 1:5)])

# ---------------------------------------------------------------------------
# test_stability
# ---------------------------------------------------------------------------

test_that("test_stability returns correct structure", {
  res <- test_stability(R, dat$Y, dat$S_e, dat$S_o)
  expect_s3_class(res, "stability_test")
  expect_s3_class(res, "data.frame")
  expect_named(res, c("y_level", "r_col", "n_exp", "n_obs", "ks_stat", "p_value"))
  # 2 y_levels × 5 R columns = 10 rows
  expect_equal(nrow(res), 10L)
})

test_that("test_stability produces finite KS stats (enough data)", {
  res <- test_stability(R, dat$Y, dat$S_e, dat$S_o)
  expect_true(all(is.finite(res$ks_stat)))
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
})

test_that("test_stability respects y_levels argument", {
  res <- test_stability(R, dat$Y, dat$S_e, dat$S_o, y_levels = 0)
  expect_equal(nrow(res), 5L)
  expect_true(all(res$y_level == 0))
})

test_that("test_stability errors when Y not observed in S_e", {
  # No overlap: all experimental Y is NA
  dat0 <- sim_rsv_data(n_e = 100, n_o = 200, n_v = 0, seed = 1)
  R0   <- as.matrix(dat0[, paste0("R", 1:5)])
  expect_error(
    test_stability(R0, dat0$Y, dat0$S_e, dat0$S_o),
    "not observed in the experimental sample"
  )
})

test_that("test_stability errors when Y not observed in S_o", {
  # Manually construct: non-overlapping design, Y in S_e but not S_o
  n <- 200
  S_e_m <- c(rep(1L, 100), rep(0L, 100))
  S_o_m <- c(rep(0L, 100), rep(1L, 100))
  Y_m   <- c(sample(0:1, 100, replace = TRUE), rep(NA_integer_, 100))
  R_m   <- matrix(rnorm(n * 5), n, 5)
  expect_error(
    test_stability(R_m, Y_m, S_e_m, S_o_m),
    "not observed in the observational sample"
  )
})

test_that("test_stability errors on dimension mismatch", {
  expect_error(test_stability(R[-1, ], dat$Y, dat$S_e, dat$S_o), "same number of rows")
  expect_error(test_stability(R, dat$Y, dat$S_e[-1], dat$S_o), "same length")
})

test_that("test_stability gives NA stats when group has < 2 obs", {
  res <- test_stability(R, dat$Y, dat$S_e, dat$S_o, y_levels = 99)
  expect_true(all(is.na(res$ks_stat)))
  expect_true(all(is.na(res$p_value)))
})

test_that("print.stability_test produces output", {
  res <- test_stability(R, dat$Y, dat$S_e, dat$S_o)
  out <- capture.output(print(res))
  expect_true(any(grepl("Stability", out)))
  expect_true(any(grepl("KS", out, ignore.case = TRUE)))
})

# ---------------------------------------------------------------------------
# test_no_direct_effect
# ---------------------------------------------------------------------------

test_that("test_no_direct_effect returns correct structure", {
  res <- test_no_direct_effect(R, dat$Y, dat$D, dat$S_e)
  expect_s3_class(res, "no_direct_effect_test")
  expect_s3_class(res, "data.frame")
  expect_named(res, c("y_level", "r_col", "n_d0", "n_d1", "ks_stat", "p_value"))
  expect_equal(nrow(res), 10L)
})

test_that("test_no_direct_effect produces finite KS stats (enough data)", {
  res <- test_no_direct_effect(R, dat$Y, dat$D, dat$S_e)
  expect_true(all(is.finite(res$ks_stat)))
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
})

test_that("test_no_direct_effect respects y_levels argument", {
  res <- test_no_direct_effect(R, dat$Y, dat$D, dat$S_e, y_levels = 1)
  expect_equal(nrow(res), 5L)
  expect_true(all(res$y_level == 1))
})

test_that("test_no_direct_effect errors when Y not observed in S_e", {
  dat0 <- sim_rsv_data(n_e = 100, n_o = 200, n_v = 0, seed = 1)
  R0   <- as.matrix(dat0[, paste0("R", 1:5)])
  expect_error(
    test_no_direct_effect(R0, dat0$Y, dat0$D, dat0$S_e),
    "not observed in the experimental sample"
  )
})

test_that("test_no_direct_effect errors when only one treatment arm", {
  D_mod <- dat$D
  D_mod[dat$S_e == 1] <- 0L
  expect_error(
    test_no_direct_effect(R, dat$Y, D_mod, dat$S_e),
    "Both D = 0 and D = 1 must be present"
  )
})

test_that("test_no_direct_effect errors on dimension mismatch", {
  expect_error(test_no_direct_effect(R[-1, ], dat$Y, dat$D, dat$S_e), "same number of rows")
  expect_error(test_no_direct_effect(R, dat$Y, dat$D[-1], dat$S_e), "same length")
})

test_that("test_no_direct_effect gives NA stats for nonexistent y_level", {
  res <- test_no_direct_effect(R, dat$Y, dat$D, dat$S_e, y_levels = 99)
  expect_true(all(is.na(res$ks_stat)))
})

test_that("print.no_direct_effect_test produces output", {
  res <- test_no_direct_effect(R, dat$Y, dat$D, dat$S_e)
  out <- capture.output(print(res))
  expect_true(any(grepl("direct.effect", out, ignore.case = TRUE)))
})
