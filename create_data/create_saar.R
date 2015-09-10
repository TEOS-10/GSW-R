library(ncdf4)
nc <- nc_open("~/src/GSW-Fortran/test/gsw_data_v3_0.nc")
p_ref <- ncvar_get(nc, "p_ref")
lats_ref <- ncvar_get(nc, "lats_ref")
longs_ref <- ncvar_get(nc, "longs_ref")
ndepth_ref <- ncvar_get(nc, "ndepth_ref")
saar_ref <- ncvar_get(nc, "SAAR_ref")
delta_sa_ref <- ncvar_get(nc, "SA_ref")
saar <- list(gsw_nx=length(longs_ref), gsw_ny=length(lats_ref), gsw_nz=length(p_ref),
             longs_ref=longs_ref, lats_ref=lats_ref, p_ref=p_ref, ndepth_ref=ndepth_ref,
             saar_ref=saar_ref, delta_sa_ref=delta_sa_ref)
save(saar, file="saar.rda")
tools::resaveRdaFiles("saar.rda")
nc_close(nc)

