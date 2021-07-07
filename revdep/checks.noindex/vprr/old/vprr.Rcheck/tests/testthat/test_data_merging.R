context("merge data frames")

library(vprr)

data('ctd_dat_combine')
data('roi_dat_combine')
data("roimeas_dat_combine")


test_that("CTD and ROI data merges",{

  expect_message(ctd_roi_merge <- vpr_ctdroi_merge(ctd_dat_combine, roi_dat_combine))
  expect_equal(length(ctd_roi_merge[[1]]), length(ctd_dat_combine[[1]])) #check that all data is merged
  expect_equal(length(ctd_roi_merge), 28) #check that all vars are present


  expect_equal(range(ctd_roi_merge$time_ms), range(ctd_dat_combine$time_ms)) #check that time range matches
  expect_true(min(roi_dat_combine$time_ms) > min(ctd_roi_merge$time_ms)) # check that all ROI data is within CTD time range
  # expect_true(max(roi_dat_combine$time_ms) < max(ctd_roi_merge$time_ms)) # not true for sample data because of subsetting

  # check for NAs
  expect_false(anyNA(ctd_roi_merge$time_ms))
  expect_false(anyNA(ctd_roi_merge$depth))
  expect_false(anyNA(ctd_roi_merge$Calanus))
  expect_false(anyNA(ctd_roi_merge$n_roi_total))


})
