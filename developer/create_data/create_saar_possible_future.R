rm(list=ls())

library(ncdf4)
nc <- nc_open("~/git/GSW-Fortran/test/gsw_data_v3_0.nc")
## Use as.vector() since these will all get handed into C, which does not understand matrices.
p_ref <- as.vector(ncvar_get(nc, "p_ref"))
lats_ref <- as.vector(ncvar_get(nc, "lats_ref"))
longs_ref <- as.vector(ncvar_get(nc, "longs_ref"))
ndepth_ref <- as.vector(ncvar_get(nc, "ndepth_ref"))
ndepth_ref[!is.finite(ndepth_ref)] <- -9e99
saar_ref <- as.vector(ncvar_get(nc, "SAAR_ref"))
saar_ref[!is.finite(saar_ref)] <- -9e99

delta_sa_ref <- as.vector(ncvar_get(nc, "deltaSA_ref"))
delta_sa_ref[!is.finite(delta_sa_ref)] <- -9e99

saar <- list(gsw_nx=length(longs_ref), gsw_ny=length(lats_ref), gsw_nz=length(p_ref),
             longs_ref=longs_ref, lats_ref=lats_ref, p_ref=p_ref, ndepth_ref=ndepth_ref,
             saar_ref=saar_ref, delta_sa_ref=delta_sa_ref)
save(saar, file="saar.rda")
tools::resaveRdaFiles("saar.rda")
nc_close(nc)

# Check the netcdf file against the matlab file
library(R.matlab)
m <- readMat("~/src/gsw_matlab_v3_06_12/library/gsw_data_v3_0.mat")
m$p_ref <- as.vector(m$p.ref)
stopifnot(all.equal(p_ref, m$p_ref))

m$lats_ref <- as.vector(m$lats.ref)
stopifnot(all.equal(lats_ref, m$lats_ref))

m$longs_ref <- as.vector(m$longs.ref)
stopifnot(all.equal(longs_ref, m$longs_ref))

m$ndepth_ref <- as.vector(m$ndepth.ref)
m$ndepth_ref[!is.finite(m$ndepth_ref)] <- -9e99
stopifnot(all.equal(ndepth_ref, m$ndepth_ref))

m$saar_ref <- as.vector(m$SAAR.ref)
m$saar_ref[!is.finite(m$saar_ref)] <- -9e99
all.equal(saar_ref, m$saar_ref) # mean relative difference 1.153597
quantile(saar_ref)
quantile(m$saar_ref)

m$delta_sa_ref <- as.vector(m$deltaSA.ref)
m$delta_sa_ref[!is.finite(m$delta_sa_ref)] <- -9e99
stopifnot(all.equal(delta_sa_ref, m$delta_sa_ref))

# So, no delta-SA differences, but we do have SAAR discrepancies, so let's
# graph them.  They are a few percent.

if (!interactive()) png("create_saar_new.png", unit="in", width=7, height=5, res=150)
par(mfrow=c(2,1), mar=c(3,3,3,1), mgp=c(2,0.7,0))
look_saar_ref <- saar_ref != -9e99
look_msaar_ref <- m$saar_ref != -9e99
look <- look_saar_ref & look_msaar_ref
saar_ref_diff <- (saar_ref - m$saar_ref)[look]
#jboxplot(saar_ref_diff)
hist(saar_ref_diff, breaks=9000, xlim=2e-5*c(-1,1), xlab="Diff.", main="")
mtext(sprintf("Histogram of saar_ref diff from NetCDF to Matlab\nNB: %.1f%% of NetCDF are missing: %.1f%% of Matlab are missing",
              100*sum(!look_saar_ref)/length(look_saar_ref),
              100*sum(!look_msaar_ref)/length(look_msaar_ref)))
hist(saar_ref[look_saar_ref], breaks=200, xlab="saar_ref", main="")

if (!interactive()) dev.off()

