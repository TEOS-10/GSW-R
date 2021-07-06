# Check the netcdf file against the matlab file
library(R.matlab)
m <- readMat("~/src/gsw_matlab_v3_06_11/library/gsw_data_v3_0.mat")
m$p_ref <- as.vector(m$p.ref)

m$lats_ref <- as.vector(m$lats.ref)

m$longs_ref <- as.vector(m$longs.ref)

m$ndepth_ref <- as.vector(m$ndepth.ref)
m$ndepth_ref[!is.finite(m$ndepth_ref)] <- -9e99

m$saar_ref <- as.vector(m$SAAR.ref)
m$saar_ref[!is.finite(m$saar_ref)] <- -9e99

m$delta_sa_ref <- as.vector(m$deltaSA.ref)
m$delta_sa_ref[!is.finite(m$delta_sa_ref)] <- -9e99

saar <- list(gsw_nx=length(m$longs_ref), gsw_ny=length(m$lats_ref), gsw_nz=length(m$p_ref),
              longs_ref=m$longs_ref, lats_ref=m$lats_ref, p_ref=m$p_ref, ndepth_ref=m$ndepth_ref,
              saar_ref=m$saar_ref, delta_sa_ref=m$delta_sa_ref)
save(saar, file="saar.rda")
tools::resaveRdaFiles("saar.rda")

