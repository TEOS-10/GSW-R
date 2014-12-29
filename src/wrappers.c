#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <R.h>
#include <Rdefines.h>
#include "gswteos-10.h"

//GSW_INVALID_VALUE	

// PART 1: macros for wrappers
//
// Any reasonable C programmer should understand how these macros work,
// in light of Part 2. The number in the macro name is the number of 
// arguments sent to the C library function named 'cname'. See part 3
// for the handling of C library functions that return void.
#define W2(wname, cname, arg1, arg2, n, rval) \
void (wname)(double *(arg1), double *(arg2), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        (rval)[i] = (cname)((arg1)[i], (arg2)[i]);\
        if ((rval)[i] == GSW_INVALID_VALUE) {\
            (rval)[i] = NA_REAL;\
        }\
    }\
}

#define W3(wname, cname, arg1, arg2, arg3, n, rval) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        (rval)[i] = (cname)((arg1)[i], (arg2)[i], (arg3)[i]);\
        if ((rval)[i] == GSW_INVALID_VALUE) {\
            (rval)[i] = NA_REAL;\
        }\
    }\
}

#define W4(wname, cname, arg1, arg2, arg3, arg4, n, rval) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), double *(arg4), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        (rval)[i] = (cname)((arg1)[i], (arg2)[i], (arg3)[i], (arg4)[i]);\
        if ((rval)[i] == GSW_INVALID_VALUE) {\
            (rval)[i] = NA_REAL;\
        }\
    }\
}

// PART 2: wrappers for functions that return a value. Wrapping is necessary 
// because the R function .C() cannot handle return values.
// See Part 3 for functios returning void.
//
// The case in the wrapper function names (first argument to macros)
// should match that in the TEOS-10 matlab functions, but note that case in
// the second argument should be lower-case, since this is used in the C
// library.
W3(wrap_gsw_adiabatic_lapse_rate_from_CT, gsw_adiabatic_lapse_rate_from_ct, SA, CT, p, n, rval)
W3(wrap_gsw_alpha, gsw_alpha, SA, CT, p, n, rval)
W3(wrap_gsw_alpha_on_beta, gsw_alpha_on_beta, SA, CT, p, n, rval)
W3(wrap_gsw_alpha_wrt_t_exact, gsw_alpha_wrt_t_exact, SA, t, p, n, rval)
W3(wrap_gsw_beta, gsw_beta, SA, CT, p, n, rval)
W3(wrap_gsw_beta_const_t_exact, gsw_beta_const_t_exact, SA, t, p, n, rval)
W3(wrap_gsw_C_from_SP, gsw_c_from_sp, SP, t, p, n, rval)
W3(wrap_gsw_cabbeling, gsw_cabbeling, SA, CT, p, n, rval)
W3(wrap_gsw_CT_freezing, gsw_ct_freezing, SA, p, saturation_fraction, n, rval)
W2(wrap_gsw_CT_from_pt, gsw_ct_from_pt, SA, pt, n, rval)
W3(wrap_gsw_cp_t_exact, gsw_cp_t_exact, SA, t, p, n, rval)
W3(wrap_gsw_CT_from_t, gsw_ct_from_t, SA, t, p, n, rval)
W4(wrap_gsw_deltaSA_from_SP, gsw_deltasa_from_sp, SP, p, longitude, latitude, n, rval)
W3(wrap_gsw_dynamic_enthalpy, gsw_dynamic_enthalpy, SA, CT, p, n, rval)
W3(wrap_gsw_enthalpy, gsw_enthalpy, SA, CT, p, n, rval)
W3(wrap_gsw_enthalpy_t_exact, gsw_enthalpy_t_exact, SA, t, p, n, rval)
// gsw_entropy_from_CT is not in the C library.
W3(wrap_gsw_entropy_from_t, gsw_entropy_from_t, SA, t, p, n, rval)
W2(wrap_gsw_grav, gsw_grav, latitude, p, n, rval)
W3(wrap_gsw_internal_energy, gsw_internal_energy, SA, CT, p, n, rval)
W2(wrap_gsw_latentheat_evap_CT, gsw_latentheat_evap_ct, SA, CT, n, rval)
W2(wrap_gsw_latentheat_evap_t, gsw_latentheat_evap_t, SA, t, n, rval)
W2(wrap_gsw_latentheat_melting, gsw_latentheat_melting, SA, p, n, rval)
W4(wrap_gsw_pot_rho_t_exact, gsw_pot_rho_t_exact, SA, t, p, p_ref, n, rval)
W3(wrap_gsw_pt0_from_t, gsw_pt0_from_t, SA, t, p, n, rval)
W2(wrap_gsw_pt_from_CT, gsw_pt_from_ct, SA, CT, n, rval)
W4(wrap_gsw_pt_from_t, gsw_pt_from_t, SA, t, p, p_ref, n, rval)
W3(wrap_gsw_rho, gsw_rho, SA, CT, p, n, rval)
W3(wrap_gsw_rho_t_exact, gsw_rho_t_exact, SA, t, p, n, rval)
W3(wrap_gsw_SA_from_rho, gsw_sa_from_rho, rho, CT, p, n, rval)
W4(wrap_gsw_SA_from_SP, gsw_sa_from_sp, CT, p, longitude, latitude, n, rval)
W4(wrap_gsw_SA_from_Sstar, gsw_sa_from_sstar, Sstar, p, longitude, latitude, n, rval)
W2(wrap_gsw_sigma0, gsw_sigma0, SA, CT, n, rval)
W2(wrap_gsw_sigma1, gsw_sigma1, SA, CT, n, rval)
W2(wrap_gsw_sigma2, gsw_sigma2, SA, CT, n, rval)
W2(wrap_gsw_sigma3, gsw_sigma3, SA, CT, n, rval)
W2(wrap_gsw_sigma4, gsw_sigma4, SA, CT, n, rval)
W3(wrap_gsw_sound_speed, gsw_sound_speed, SA, t, p, n, rval)
W3(wrap_gsw_sound_speed_t_exact, gsw_sound_speed_t_exact, SA, t, p, n, rval)
// gsw_specvol coded in R
W3(wrap_gsw_specvol_anom, gsw_specvol_anom, SA, CT, p, n, rval)
W3(wrap_gsw_specvol_t_exact, gsw_specvol_t_exact, SA, t, p, n, rval)
W3(wrap_gsw_SP_from_C, gsw_sp_from_c, C, t, p, n, rval)
W4(wrap_gsw_SP_from_SA, gsw_sp_from_sa, SA, p, longitude, latitude, n, rval)
W3(wrap_gsw_t_freezing, gsw_t_freezing, SA, p, saturation_fraction, n, rval)
W3(wrap_gsw_t_from_CT, gsw_t_from_ct, SA, CT, p, n, rval)
W3(wrap_gsw_thermobaric, gsw_thermobaric, SA, CT, p, n, rval)
W2(wrap_gsw_z_from_p, gsw_z_from_p, p, lat, n, rval)

// PART 3
//
// Functions returning void do not really need wrappers, but the R code
// is simpler to read if we use them anyway.
void wrap_gsw_rho_first_derivatives(double *SA, double *CT, double *p, int *n,
        double *drho_dsa, double *drho_dct, double *drho_dp)
{
    extern void gsw_rho_first_derivatives(double sa, double ct, double p,
            double *drho_dsa, double *drho_dct, double *drho_dp);
    for (int i=0; i < *(n); i++) {
        gsw_rho_first_derivatives(SA[i], CT[i], p[i], &drho_dsa[i], &drho_dct[i], &drho_dp[i]);
    }
}

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

