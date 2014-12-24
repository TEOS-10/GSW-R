#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "gswteos-10.h"

// PART 1: macros for wrappers
//
// Any reasonable C programmer should understand how these macros work,
// in light of Part 2. 
//
// wrap a 2-arg function that returns 1 value
#define W21(wname, cname, arg1, arg2, n, rval) \
void (wname)(double *(arg1), double *(arg2), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        (rval)[i] = (cname)((arg1)[i], (arg2)[i]);\
    }\
}
//
// wrap a 3-arg function that returns 1 value
#define W31(wname, cname, arg1, arg2, arg3, n, rval) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        (rval)[i] = (cname)((arg1)[i], (arg2)[i], (arg3)[i]);\
    }\
}
//
// wrap a 4-arg function that returns 1 value
#define W41(wname, cname, arg1, arg2, arg3, arg4, n, rval) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), double *(arg4), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        (rval)[i] = (cname)((arg1)[i], (arg2)[i], (arg3)[i], (arg4)[i]);\
    }\
}

// PART 2: wrappers for functions that return a value (which cannot be called
// directly with .C() in R). See Part 3 for functios returning void.
//
// Note
// that the case in the wrapper function names (first argument to macros)
// should match that in the TEOS-10 matlab functions, but note that case in
// the second argument should be lower-case, since this is used in the C
// library.
//
// Please put functions into alphabetical order here, ignoring case.
W31(wrap_gsw_adiabatic_lapse_rate_from_CT, gsw_adiabatic_lapse_rate_from_ct, SA, CT, p, n, rval)
W31(wrap_gsw_alpha, gsw_alpha, SA, CT, p, n, rval)
W31(wrap_gsw_alpha_on_beta, gsw_alpha_on_beta, SA, CT, p, n, rval)
W31(wrap_gsw_alpha_wrt_t_exact, gsw_alpha_wrt_t_exact, SA, t, p, n, rval)
W31(wrap_gsw_beta, gsw_beta, SA, CT, p, n, rval)
W31(wrap_gsw_C_from_SP, gsw_c_from_sp, SP, t, p, n, rval)
W31(wrap_gsw_cabbeling, gsw_cabbeling, SA, CT, p, n, rval)
W31(wrap_gsw_CT_freezing, gsw_ct_freezing, SA, p, saturation_fraction, n, rval)
W21(wrap_gsw_CT_from_pt, gsw_ct_from_pt, SA, pt, n, rval)
W31(wrap_gsw_cp_t_exact, gsw_cp_t_exact, SA, t, p, n, rval)
W31(wrap_gsw_CT_from_t, gsw_ct_from_t, SA, t, p, n, rval)
W21(wrap_gsw_grav, gsw_grav, latitude, p, n, rval)
W41(wrap_gsw_pot_rho_t_exact, gsw_pot_rho_t_exact, SA, t, p, p_ref, n, rval)
W31(wrap_gsw_rho, gsw_rho, SA, CT, p, n, rval)
W31(wrap_gsw_rho_t_exact, gsw_rho_t_exact, SA, t, p, n, rval)
W41(wrap_gsw_SA_from_SP, gsw_sa_from_sp, SA, p, longitude, latitude, n, rval)
W21(wrap_gsw_sigma0, gsw_sigma0, SA, CT, n, rval)
W21(wrap_gsw_sigma1, gsw_sigma1, SA, CT, n, rval)
W21(wrap_gsw_sigma2, gsw_sigma2, SA, CT, n, rval)
W21(wrap_gsw_sigma3, gsw_sigma3, SA, CT, n, rval)
W21(wrap_gsw_sigma4, gsw_sigma4, SA, CT, n, rval)
W31(wrap_gsw_sound_speed, gsw_sound_speed, SA, t, p, n, rval)
W31(wrap_gsw_sound_speed_t_exact, gsw_sound_speed_t_exact, SA, t, p, n, rval)
W31(wrap_gsw_SP_from_C, gsw_sp_from_c, C, t, p, n, rval)
W41(wrap_gsw_SP_from_SA, gsw_sp_from_sa, SA, p, longitude, latitude, n, rval)
W31(wrap_gsw_t_freezing, gsw_t_freezing, SA, p, saturation_fraction, n, rval)
W31(wrap_gsw_t_from_CT, gsw_t_from_ct, SA, CT, p, n, rval)
W21(wrap_gsw_z_from_p, gsw_z_from_p, p, latitude, n, rval)

// PART 3
//
// Functions returning void do not really need wrappers, but the R code
// is simpler to read if we use them anyway.
void wrap_gsw_Nsquared(double *SA, double *CT, double *p, double *latitude, int *n, double *n2, double *p_mid)
{
    extern void gsw_nsquared(double *sa, double *ct, double *p, double *latitude, int nz, double *n2, double *p_mid);
    gsw_nsquared(SA, CT, p, latitude, *n, n2, p_mid);
}

void wrap_gsw_Turner_Rsubrho(double *SA, double *CT, double *p, int *n, double *Tu, double *Rsubrho, double *p_mid)
{
    extern void gsw_turner_rsubrho(double *sa, double *ct, double *p, int nz, double *Tu, double *Rsubrho, double *p_mid);
    gsw_turner_rsubrho(SA, CT, p, *n, Tu, Rsubrho, p_mid);
}

