#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "gswteos-10.h"

// Case in wrapper function names should match that in the TEOS-10 matlab functions,
// but note that the gsw_ functions we are calling all have lower-case names.
// Please put functions into alphabetical order here, ignoring case.

void wrap_gsw_CT_from_t(double *SA, double *t, double *p, int *n, double *rval)
{
    extern double gsw_ct_from_t(double SA, double t, double p);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_ct_from_t(SA[i], t[i], p[i]);
}

void wrap_gsw_Nsquared(double *SA, double *CT, double *p, double *latitude, int *n, double *n2, double *p_mid)
{
    extern void gsw_nsquared(double *sa, double *ct, double *p, double *latitude, int nz, double *n2, double *p_mid);
    gsw_nsquared(SA, CT, p, latitude, *n, n2, p_mid);
}

void wrap_gsw_SP_from_C(double *C, double *t, double *p, int *n, double *rval)
{
    extern double gsw_SP_from_c(double C, double t, double p);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_sp_from_c(C[i], t[i], p[i]);
}

void wrap_gsw_SA_from_SP(double *SA, double *p, double *longitude, double *latitude,
        int *n, double *rval)
{
    extern double gsw_SA_from_SP(double SA, double p, double longitude, double latitude);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_sa_from_sp(SA[i], p[i], longitude[i], latitude[i]);
}


