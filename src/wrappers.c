#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "gswteos-10.h"

// Case in wrapper function names should match that in the TEOS-10 matlab functions,
// but note that the gsw_ functions we are calling all have lower-case names.
// Please put functions into alphabetical order here, ignoring case.

void wrap_gsw_alpha(double *SA, double *CT, double *p, int *n, double *rval)
{
    extern double gsw_alpha(double SA, double CT, double p);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_alpha(SA[i], CT[i], p[i]);
}

void wrap_gsw_beta(double *SA, double *CT, double *p, int *n, double *rval)
{
    extern double gsw_beta(double SA, double CT, double p);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_beta(SA[i], CT[i], p[i]);
}

void wrap_gsw_C_from_SP(double *SP, double *t, double *p, int *n, double *rval)
{
    extern double gsw_c_from_sp(double SP, double t, double p);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_c_from_sp(SP[i], t[i], p[i]);
}

void wrap_gsw_CT_from_t(double *SA, double *t, double *p, int *n, double *rval)
{
    extern double gsw_ct_from_t(double SA, double t, double p);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_ct_from_t(SA[i], t[i], p[i]);
}

void wrap_gsw_grav(double *latitude, double *p, int *n, double *rval)
{
    extern double gsw_grav(double latitude, double p);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_grav(latitude[i], p[i]);
}

void wrap_gsw_Nsquared(double *SA, double *CT, double *p, double *latitude, int *n, double *n2, double *p_mid)
{
    extern void gsw_nsquared(double *sa, double *ct, double *p, double *latitude, int nz, double *n2, double *p_mid);
    gsw_nsquared(SA, CT, p, latitude, *n, n2, p_mid);
}

void wrap_gsw_pot_rho_t_exact(double *SA, double *t, double *p, double *p_ref, int *n, double *rval)
{
    extern double gsw_pot_rho_t_exact(double sa, double t, double p, double p_ref);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_pot_rho_t_exact(SA[i], t[i], p[i], p_ref[i]);
}

void wrap_gsw_rho(double *SA, double *CT, double *p, int *n, double *rval)
{
    extern double gsw_rho(double SA, double CT, double p);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_rho(SA[i], CT[i], p[i]);
}

void wrap_gsw_rho_t_exact(double *SA, double *t, double *p, int *n, double *rval)
{
    extern double gsw_rho_t_exact(double SA, double t, double p);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_rho_t_exact(SA[i], t[i], p[i]);
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

void wrap_gsw_sigma0(double *SA, double *CT,
        int *n, double *rval)
{
    extern double gsw_sigma0(double SA, double CT);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_sigma0(SA[i], CT[i]);
}

void wrap_gsw_sigma1(double *SA, double *CT,
        int *n, double *rval)
{
    extern double gsw_sigma1(double SA, double CT);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_sigma1(SA[i], CT[i]);
}

void wrap_gsw_sigma2(double *SA, double *CT,
        int *n, double *rval)
{
    extern double gsw_sigma2(double SA, double CT);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_sigma2(SA[i], CT[i]);
}

void wrap_gsw_sigma3(double *SA, double *CT,
        int *n, double *rval)
{
    extern double gsw_sigma3(double SA, double CT);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_sigma3(SA[i], CT[i]);
}

void wrap_gsw_sigma4(double *SA, double *CT,
        int *n, double *rval)
{
    extern double gsw_sigma4(double SA, double CT);
    for (int i=0; i < *n; i++)
        rval[i] = gsw_sigma4(SA[i], CT[i]);
}

void wrap_gsw_Turner_Rsubrho(double *SA, double *CT, double *p, int *n, double *Tu, double *Rsubrho, double *p_mid)
{
    extern void gsw_turner_rsubrho(double *sa, double *ct, double *p, int nz, double *Tu, double *Rsubrho, double *p_mid);
    gsw_turner_rsubrho(SA, CT, p, *n, Tu, Rsubrho, p_mid);
}

void wrap_gsw_z_from_p(double *p, double *latitude, int *n, double *rval)
{
    extern double gsw_z_from_p(double p, double latitude);
    for (int i=0; i < *n; i++)
      rval[i] = gsw_z_from_p(p[i], latitude[i]);
}

