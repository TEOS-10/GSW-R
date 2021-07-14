pkgname <- "vprr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('vprr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("vpr_autoid_read")
### * vpr_autoid_read

flush(stderr()); flush(stdout())

### Name: vpr_autoid_read
### Title: Read VPR aid files
### Aliases: vpr_autoid_read

### ** Examples


station_of_interest <- 'test'
dayhour <- c('d222.h03', 'd222.h04')

#' #VPR OPTICAL SETTING (S0, S1, S2 OR S3)
opticalSetting <- "S2"
imageVolume <- 83663 #mm^3

auto_id_folder <- system.file('extdata/COR2019002/autoid/', package = 'vprr', mustWork = TRUE)
auto_id_path <- list.files(paste0(auto_id_folder, "/"), full.names = TRUE)

#'   # Path to aid for each taxa
aid_path <- paste0(auto_id_path, '/aid/')
# Path to mea for each taxa
aidmea_path <- paste0(auto_id_path, '/aidmea/')

# AUTO ID FILES
aid_file_list <- list()
aidmea_file_list <- list()
for (i in 1:length(dayhour)) {
  aid_file_list[[i]] <-
    list.files(aid_path, pattern = dayhour[[i]], full.names = TRUE)
  # SIZE DATA FILES
  aidmea_file_list[[i]] <-
    list.files(aidmea_path, pattern = dayhour[[i]], full.names = TRUE)
}

aid_file_list_all <- unlist(aid_file_list)
aidmea_file_list_all <- unlist(aidmea_file_list)

 # ROIs
roi_dat_combine <-
  vpr_autoid_read(
    file_list_aid = aid_file_list_all,
    file_list_aidmeas = aidmea_file_list_all,
    export = 'aid',
    station_of_interest = station_of_interest,
    opticalSetting = opticalSetting,
    warn = FALSE
  )

# MEASUREMENTS
roimeas_dat_combine <-
  vpr_autoid_read(
    file_list_aid = aid_file_list_all,
    file_list_aidmeas = aidmea_file_list_all,
    export = 'aidmeas',
    station_of_interest = station_of_interest,
    opticalSetting = opticalSetting,
    warn = FALSE
 )




cleanEx()
nameEx("vpr_category")
### * vpr_category

flush(stderr()); flush(stdout())

### Name: vpr_category
### Title: Get taxa ids from string
### Aliases: vpr_category

### ** Examples

taxa_string <- 'C:/data/cruise/autoid/Calanus/d000/h00'
vpr_category(taxa_string)




cleanEx()
nameEx("vpr_ctd_read")
### * vpr_ctd_read

flush(stderr()); flush(stdout())

### Name: vpr_ctd_read
### Title: Read and format CTD VPR data
### Aliases: vpr_ctd_read

### ** Examples


station_of_interest <- 'test'

ctd_files <- system.file("extdata/COR2019002/rois/vpr5/d222", "h03ctd.dat",
package = "vprr", mustWork = TRUE)

ctd_dat_combine <- vpr_ctd_read(ctd_files, station_of_interest)




cleanEx()
nameEx("vpr_ctd_ymd")
### * vpr_ctd_ymd

flush(stderr()); flush(stdout())

### Name: vpr_ctd_ymd
### Title: Add Year/ month/ day hour:minute:second information
### Aliases: vpr_ctd_ymd

### ** Examples

year <- 2019
data('ctd_roi_merge')
dat <- vpr_ctd_ymd(ctd_roi_merge, year)





cleanEx()
nameEx("vpr_ctdroi_merge")
### * vpr_ctdroi_merge

flush(stderr()); flush(stdout())

### Name: vpr_ctdroi_merge
### Title: Merge CTD and ROI data from VPR
### Aliases: vpr_ctdroi_merge

### ** Examples

data('ctd_dat_combine')
data('roi_dat_combine')

ctd_roi_merge <- vpr_ctdroi_merge(ctd_dat_combine, roi_dat_combine)



cleanEx()
nameEx("vpr_ctdroisize_merge")
### * vpr_ctdroisize_merge

flush(stderr()); flush(stdout())

### Name: vpr_ctdroisize_merge
### Title: Format CTD and Size data from VPR
### Aliases: vpr_ctdroisize_merge

### ** Examples


data("ctd_roi_merge")
data("roimeas_dat_combine")
category_of_interest = 'Calanus'

ctd_roi_merge$avg_hr <- ctd_roi_merge$time_ms /3.6e+06

size_df_f <- vpr_ctdroisize_merge(ctd_roi_merge, data_mea = roimeas_dat_combine,
 taxa_of_interest = category_of_interest)




cleanEx()
nameEx("vpr_day")
### * vpr_day

flush(stderr()); flush(stdout())

### Name: vpr_day
### Title: Get day identifier
### Aliases: vpr_day

### ** Examples

day_string <- 'C:/data/cruise/autoid/Calanus/d000/h00'
vpr_day(day_string)




cleanEx()
nameEx("vpr_hour")
### * vpr_hour

flush(stderr()); flush(stdout())

### Name: vpr_hour
### Title: Get hour identifier
### Aliases: vpr_hour

### ** Examples

hour_string <- 'C:/data/cruise/autoid/Calanus/d000/h00'
vpr_hour(hour_string)




cleanEx()
nameEx("vpr_oce_create")
### * vpr_oce_create

flush(stderr()); flush(stdout())

### Name: vpr_oce_create
### Title: Create ctd oce object with vpr data
### Aliases: vpr_oce_create

### ** Examples

data('ctd_roi_merge')
oce_dat <- vpr_oce_create(ctd_roi_merge)




cleanEx()
nameEx("vpr_roi")
### * vpr_roi

flush(stderr()); flush(stdout())

### Name: vpr_roi
### Title: Get roi ids from string
### Aliases: vpr_roi

### ** Examples


roi_string <- 'roi.0100000000.tif'
vpr_roi(roi_string)




cleanEx()
nameEx("vpr_roi_concentration")
### * vpr_roi_concentration

flush(stderr()); flush(stdout())

### Name: vpr_roi_concentration
### Title: Calculate VPR concentrations
### Aliases: vpr_roi_concentration

### ** Examples


data('ctd_roi_merge')
ctd_roi_merge$avg_hr <- ctd_roi_merge$time_ms /3.6e+06

taxas_list <- c('Calanus', 'krill')
binSize <- 5
station_of_interest <- 'test'
imageVolume <- 83663

taxa_conc_n <- vpr_roi_concentration(ctd_roi_merge, taxas_list,
station_of_interest, binSize, imageVolume)




cleanEx()
nameEx("vpr_save")
### * vpr_save

flush(stderr()); flush(stdout())

### Name: vpr_save
### Title: Save VPR data as an as.oce object
### Aliases: vpr_save

### ** Examples

data("taxa_conc_n")
metadata <- c('deploymentType' = 'towyo', 'waterDepth' =
max(ctd_roi_merge$pressure), 'serialNumber' = NA, 'latitude' = 47,
'longitude' = -65, 'castDate' = '2019-08-11', 'castStartTime'= '00:00',
'castEndTime' = '01:00', 'processedBy' = 'E. Chisholm', 'opticalSetting' =
'S2', 'imageVolume' = 83663, 'comment' = 'test data')

oce_dat <- vpr_save(taxa_conc_n, metadata)
# save(oce_dat, file = vpr_save.RData') # save data




cleanEx()
nameEx("vpr_size_bin")
### * vpr_size_bin

flush(stderr()); flush(stdout())

### Name: vpr_size_bin
### Title: Bin VPR size data
### Aliases: vpr_size_bin

### ** Examples


data('size_df_f')
vpr_size_bin(size_df_f, bin_mea = 5)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
