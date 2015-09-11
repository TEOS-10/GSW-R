## NOTE: this file is for use by developers only. It is not incorporated into
## the R package, although its output ("../data/saar.rda") is.

library(ncdf4)
nc <- nc_open("~/src/GSW-Fortran/test/gsw_data_v3_0.nc")
## Use as.vector() since these will all get handed into C, which does not understand matrices.
p_ref <- as.vector(ncvar_get(nc, "p_ref"))
lats_ref <- as.vector(ncvar_get(nc, "lats_ref"))
longs_ref <- as.vector(ncvar_get(nc, "longs_ref"))
ndepth_ref <- as.vector(ncvar_get(nc, "ndepth_ref"))
ndepth_ref[!is.finite(ndepth_ref)] <- NA
saar_ref <- as.vector(ncvar_get(nc, "SAAR_ref"))
saar_ref[!is.finite(saar_ref)] <- NA
delta_sa_ref <- as.vector(ncvar_get(nc, "deltaSA_ref"))
delta_sa_ref[!is.finite(delta_sa_ref)] <- NA
saar <- list(gsw_nx=length(longs_ref), gsw_ny=length(lats_ref), gsw_nz=length(p_ref),
             longs_ref=longs_ref, lats_ref=lats_ref, p_ref=p_ref, ndepth_ref=ndepth_ref,
             saar_ref=saar_ref, delta_sa_ref=delta_sa_ref)
save(saar, file="saar.rda")
tools::resaveRdaFiles("saar.rda")
nc_close(nc)

new <- saar
load("~/src/GSW-R/data/saar.rda")
old <- saar
str(new)
str(old)

library(testthat)
expect_equal(new$gsw_nx, old$gsw_nx)
expect_equal(new$gsw_ny, old$gsw_ny)
expect_equal(new$gsw_nz, old$gsw_nz)
expect_equal(new$longs_ref, old$longs_ref)
expect_equal(new$lats_ref, old$lats_ref)
expect_equal(new$p_ref, old$p_ref)
missing <- -9e99

look <- old$ndepth_ref != missing
expect_equal(old$ndepth_ref[look], new$ndepth_ref[look])
look <- old$saar_ref != missing
expect_equal(old$saar_ref[look], new$saar_ref[look])
look <- old$delta_sa_ref != missing
expect_equal(old$delta_sa_ref[look], new$delta_sa_ref[look])

