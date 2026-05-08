test_that("sim_rsv_data returns a valid data frame", {
  dat <- sim_rsv_data(n_e = 100, n_o = 200, seed = 1)

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 300)
  expect_named(dat, c("Y", "D", "S_e", "S_o", "R1", "R2", "R3", "R4", "R5"))
})

test_that("sim_rsv_data sample indicators are correct with no overlap", {
  dat <- sim_rsv_data(n_e = 100, n_o = 200, n_v = 0, seed = 2)

  expect_equal(sum(dat$S_e), 100)
  expect_equal(sum(dat$S_o), 200)
  # No overlap: S_e and S_o are non-overlapping
  expect_equal(sum(dat$S_e == 1 & dat$S_o == 1), 0)

  # D observed iff S_e = 1
  expect_true(all(is.na(dat$D[dat$S_e == 0])))
  expect_true(all(!is.na(dat$D[dat$S_e == 1])))

  # Y in {0, 1} where observed
  expect_true(all(na.omit(dat$Y) %in% c(0L, 1L)))
})

test_that("sim_rsv_data default has Y = NA for all experimental units", {
  dat <- sim_rsv_data(n_e = 100, n_o = 200, n_v = 0, seed = 3)

  # Pure experimental: Y is NA
  expect_true(all(is.na(dat$Y[dat$S_e == 1 & dat$S_o == 0])))
  # Y observed only in observational sample
  expect_equal(sum(!is.na(dat$Y)), 200L)
})

test_that("sim_rsv_data n_v creates overlap subsample (S_e=1, S_o=1)", {
  dat <- sim_rsv_data(n_e = 100, n_o = 200, n_v = 30, seed = 4)

  # 30 overlap units
  expect_equal(sum(dat$S_e == 1 & dat$S_o == 1), 30L)
  # 70 pure experimental, 200 pure observational
  expect_equal(sum(dat$S_e == 1 & dat$S_o == 0), 70L)
  expect_equal(sum(dat$S_e == 0 & dat$S_o == 1), 200L)

  # Overlap units: Y and D both observed
  ovl <- dat$S_e == 1 & dat$S_o == 1
  expect_true(all(!is.na(dat$Y[ovl])))
  expect_true(all(!is.na(dat$D[ovl])))

  # Pure experimental: Y = NA
  pex <- dat$S_e == 1 & dat$S_o == 0
  expect_true(all(is.na(dat$Y[pex])))

  # Total observed Y = n_o + n_v
  expect_equal(sum(!is.na(dat$Y)), 200L + 30L)
})

test_that("sim_rsv_data RSVs are always observed", {
  dat <- sim_rsv_data(n_e = 100, n_o = 200, n_v = 0, seed = 5)
  r_cols <- grep("^R", names(dat), value = TRUE)
  expect_true(all(!is.na(dat[, r_cols])))
})

test_that("sim_rsv_data is reproducible with seed", {
  d1 <- sim_rsv_data(seed = 99)
  d2 <- sim_rsv_data(seed = 99)
  expect_identical(d1, d2)
})

test_that("sim_rsv_data n_r controls RSV columns", {
  dat <- sim_rsv_data(n_r = 3, seed = 7)
  expect_named(dat, c("Y", "D", "S_e", "S_o", "R1", "R2", "R3"))
})

test_that("sim_rsv_data n_v = n_e gives all experimental Y observed", {
  dat <- sim_rsv_data(n_e = 50, n_o = 100, n_v = 50, seed = 8)
  # All experimental units are overlap units
  expect_equal(sum(dat$S_e == 1 & dat$S_o == 1), 50L)
  expect_true(all(!is.na(dat$Y[dat$S_e == 1])))
})

test_that("sim_rsv_data errors on invalid n_v", {
  expect_error(sim_rsv_data(n_e = 100, n_v = 101))
  expect_error(sim_rsv_data(n_e = 100, n_v = -1))
})
