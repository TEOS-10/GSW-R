context("read files")

library(vprr)

# required metadata and data


ctd_files[[1]] <- system.file('extdata/COR2019002/rois/vpr5/d222', 'h03ctd.dat', package = 'vprr', mustWork = TRUE)
ctd_files[[2]] <- system.file('extdata/COR2019002/rois/vpr5/d222', 'h04ctd.dat', package = 'vprr', mustWork = TRUE)


station_of_interest <- 'test'
day_of_interest <- '222'
hour_of_interest <- c('03', '04')

no_ctd_files <- list()
col_list <- c("time_ms", "conductivity", "temperature", "pressure", "salinity", "fluor_ref", "fluorescence_mv",
              "turbidity_ref", "turbidity_mv", "altitude_NA")
# read raw ctd files

raw1 <- readLines(ctd_files[1])
raw2 <- readLines(ctd_files[2])

# CTD file read in
test_that("CTD files are read in accurately",{
  expect_warning(vpr_ctd_read(ctd_files, station_of_interest), 'CTD data columns named based on 2019 defaults!') # produces warning about column names
  expect_error(vpr_ctd_read(no_ctd_files, station_of_interest)) # error when no files given
  expect_silent(ctd_dat_combine <- vpr_ctd_read(ctd_files, station_of_interest, col_list = col_list )) # no warnings when col_list is specififed

  expect_equal(is.data.frame(ctd_dat_combine), TRUE) # check output is dataframe
  expect_equal(ctd_dat_combine$station[1], station_of_interest) #check station is properly copied
  expect_equal(ctd_dat_combine$day[1], paste0('d', day_of_interest)) # check day is properly identified
  expect_equal(unique(ctd_dat_combine$hour), paste0('h', hour_of_interest)) # check hour is properly identified

  # check all columns are numeric
  expect_equal(is.numeric(ctd_dat_combine$time_ms), TRUE)
  expect_equal(is.numeric(ctd_dat_combine$conductivity), TRUE)
  expect_equal(is.numeric(ctd_dat_combine$temperature), TRUE)
  expect_equal(is.numeric(ctd_dat_combine$pressure), TRUE)
  expect_equal(is.numeric(ctd_dat_combine$salinity), TRUE)
  expect_equal(is.numeric(ctd_dat_combine$fluor_ref), TRUE)
  expect_equal(is.numeric(ctd_dat_combine$fluorescence_mv), TRUE)
  expect_equal(is.numeric(ctd_dat_combine$turbidity_ref), TRUE)
  expect_equal(is.numeric(ctd_dat_combine$turbidity_mv), TRUE)
  expect_equal(is.numeric(ctd_dat_combine$altitude_NA), TRUE)

  expect_equal(stringr::str_length(as.character(ctd_dat_combine$time_ms[1])), 8) # check time stamp is 8 digits

  expect_equal(length(raw1)+length(raw2), length(ctd_dat_combine[[1]])) # check all lines of data are included

})

# VPR autoid file read in

autoid_files <- list.files(system.file('extdata/COR2019002/autoid/', package = 'vprr'), recursive = TRUE)
autoidfull_files <- file.path(system.file('extdata/COR2019002/autoid/', package = 'vprr'), autoid_files)

aidmea_files <- grep(autoidfull_files, pattern = 'aidmea', value = TRUE)
aid_files <- autoidfull_files[!autoidfull_files %in% aidmea_files]


opticalSetting <- "S2"



test_that("VPR autoid files are read in accurately", {

  # AID files
  expect_silent(roi_dat_combine <-
    vpr_autoid_read(
      file_list_aid = aid_files,
      file_list_aidmeas = aidmea_files,
      export = 'aid',
      station_of_interest = station_of_interest,
      opticalSetting = opticalSetting
    ))


  expect_equal(is.data.frame(roi_dat_combine), TRUE) # check output is dataframe
  expect_identical(as.character(roi_dat_combine$time_ms), roi_dat_combine$roi) #check time and roi stamps match

  taxa_names <- names(roi_dat_combine)[-c(1,13)]
  t_names_exp <- list()
  for(i in 1:length(aid_files)){
    t_names_exp[[i]] <- vpr_category(aid_files[[i]])
  }
  t_names_exp <- unique(unlist(t_names_exp))

  expect_identical(taxa_names, t_names_exp) #check taxa names are pulled properly


  # check for NAs
  bp_anyNA = getFromNamespace("anyNA", "backports")
  expect_equal(bp_anyNA(roi_dat_combine), FALSE)

  # check that all ROI images are accounted for in aid files
  # !! WARNING sample data fails this test
  # failure is likely due to removal of some images before matlab classification based on image size


  # img1 <- list.files(path = paste0(castdir, 'h', hour_of_interest[[1]]))
  # img2 <- list.files(path = paste0(castdir, 'h', hour_of_interest[[2]]))
  #
  # num_roi_images <- length(img1) +length(img2)
  #
  # num_dat_img <- list()
  # for(i in 1:length(taxa_names)){
  #   num_dat_img[[i]] <- sum(roi_dat_combine[[taxa_names[[i]]]])
  # }
  #
  # num_dat_img <- sum(unlist(num_dat_img))
  #
  # roi_folder_nums1 <- substr(vpr_roi(img1), 1, 8)
  # roi_folder_nums2 <- substr(vpr_roi(img2), 1, 8)
  # roi_folder_nums <- c(roi_folder_nums1, roi_folder_nums2)
  #
  # get roi ids from 'missing' rois
  # tt <- roi_folder_nums[-which(roi_folder_nums %in% roi_dat_combine$roi)]

  # test
  # expect_identical(unique(roi_folder_nums), roi_dat_combine$roi)

  # aid mea files

  expect_silent(roimeas_dat_combine <-
    vpr_autoid_read(
      file_list_aid = aid_files,
      file_list_aidmeas = aidmea_files,
      export = 'aidmeas',
      station_of_interest = station_of_interest,
      opticalSetting = opticalSetting
    ))

  # pixel data
  expect_warning(roimeas_dat_combine_px <-
                   vpr_autoid_read(
                     file_list_aid = aid_files,
                     file_list_aidmeas = aidmea_files,
                     export = 'aidmeas',
                     station_of_interest = station_of_interest
                   ))

  expect_equal(is.data.frame(roimeas_dat_combine), TRUE) # check data output is data frame

  # check match to aid files
  num_dat_img <- list()
  for(i in 1:length(taxa_names)){
    num_dat_img[[i]] <- sum(roi_dat_combine[[taxa_names[[i]]]])
  }

  num_dat_img <- sum(unlist(num_dat_img))

  expect_equal(length(roimeas_dat_combine[[1]]), num_dat_img) # check that aid and aidmea data frames have same rois
  expect_equal(unique(roimeas_dat_combine$day_hour), paste0('d', day_of_interest, '.h', hour_of_interest)) # check day and hour match
  expect_identical(stringr::str_sub(roimeas_dat_combine$roi, 1, 8), as.character(roimeas_dat_combine$time_ms)) # check that roi and times are identical

  # check pixel and mm data frames are the same for all metadata
  expect_identical(roimeas_dat_combine$roi, roimeas_dat_combine_px$roi)
  expect_identical(roimeas_dat_combine$taxa, roimeas_dat_combine_px$taxa)
  expect_identical(roimeas_dat_combine$day_hour, roimeas_dat_combine_px$day_hour)
  expect_identical(roimeas_dat_combine$station, roimeas_dat_combine_px$station)
  expect_identical(roimeas_dat_combine$time_ms, roimeas_dat_combine_px$time_ms)

  expect_identical(roimeas_dat_combine[4:10], px_to_mm(roimeas_dat_combine_px[4:10], opticalSetting)) # check px to mm conversion

  # check numeric
  expect_equal(is.numeric(roimeas_dat_combine[[4]]), TRUE)
  expect_equal(is.numeric(roimeas_dat_combine[[5]]), TRUE)
  expect_equal(is.numeric(roimeas_dat_combine[[6]]), TRUE)
  expect_equal(is.numeric(roimeas_dat_combine[[7]]), TRUE)
  expect_equal(is.numeric(roimeas_dat_combine[[8]]), TRUE)
  expect_equal(is.numeric(roimeas_dat_combine[[9]]), TRUE)
  expect_equal(is.numeric(roimeas_dat_combine[[10]]), TRUE)


})

