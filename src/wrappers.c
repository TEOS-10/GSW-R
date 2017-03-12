#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <R.h>
#include <Rdefines.h>
#include "gswteos-10.h"

// PART 1
//
// Interface to store data on package load and clear it on unload (see
// functions .onLoad() and .onUnload() in the file ../R/zzz.R).
void set_up_gsw_data(int *gsw_nx_val,
        int *gsw_ny_val,
        int *gsw_nz_val,
        double *longs_ref_val,
        double *lats_ref_val,
        double *p_ref_val,
        double *ndepth_ref_val,
        double *saar_ref_val,
        double *delta_sa_ref_val)
{
    extern int gsw_nx, gsw_ny, gsw_nz;
    extern double *longs_ref, *lats_ref, *p_ref, *ndepth_ref, *saar_ref, *delta_sa_ref;
    if (!gsw_nx) {
        //Rprintf("about to set up sarr globals\n");
        gsw_nx = *gsw_nx_val;
        gsw_ny = *gsw_ny_val;
        gsw_nz = *gsw_nz_val;
        if (!gsw_nx) error("something is wrong with gsw_nx\n");
        if (!gsw_ny) error("something is wrong with gsw_ny\n");
        if (!gsw_nz) error("something is wrong with gsw_nz\n");

        longs_ref = Calloc(gsw_nx, double);
        if (!longs_ref) error("cannot allocate memory for GSW internal data item \"longs_ref\"\n");
        for (int i = 0; i < gsw_nx; i++) longs_ref[i] = longs_ref_val[i];;

        lats_ref = Calloc(gsw_ny, double);
        if (!lats_ref) error("cannot allocate memory for GSW internal data item \"lats_ref\"\n");
        for (int i = 0; i < gsw_ny; i++) lats_ref[i] = lats_ref_val[i];;

        p_ref = Calloc(gsw_nz, double);
        if (!p_ref) error("cannot allocate memory for GSW internal data item \"p_ref\"\n");
        for (int i = 0; i < gsw_nz; i++) p_ref[i] = p_ref_val[i];;

        //  ndepth_ref isn't actually used anywhere, so perhaps it should not be saved.
        int nxy = gsw_nx * gsw_ny;
        ndepth_ref = Calloc(nxy, double);
        if (!ndepth_ref) error("cannot allocate memory for GSW internal data item \"ndepth_ref\"\n");
        for (int i = 0; i < nxy; i++) ndepth_ref[i] = ndepth_ref_val[i];;

        int nxyz = gsw_nx * gsw_ny * gsw_nz;
        saar_ref = Calloc(nxyz, double);
        if (!saar_ref) error("cannot allocate memory for GSW internal data item \"saar_ref\"\n");
        for (int i = 0; i < nxyz; i++) saar_ref[i] = saar_ref_val[i];;

        delta_sa_ref = Calloc(nxyz, double);
        if (!delta_sa_ref) error("cannot allocate memory for GSW internal data item \"delta_sa_ref\"\n");
        for (int i = 0; i < nxyz; i++) delta_sa_ref[i] = delta_sa_ref_val[i];;

        //Rprintf(" ... done setting up globals\n");
    } else {
        Rprintf("sarr globals -- already set up\n");
    }
}

void clear_gsw_data()
{
    extern int gsw_nx, gsw_ny, gsw_nz;
    extern double *longs_ref, *lats_ref, *p_ref, *ndepth_ref, *saar_ref, *delta_sa_ref;
    gsw_nx = 0;
    gsw_ny = 0;
    gsw_nz = 0;
    Free(longs_ref);
    longs_ref = NULL;
    Free(lats_ref);
    lats_ref = NULL;
    Free(p_ref);
    p_ref = NULL;
    Free(ndepth_ref);
    ndepth_ref = NULL;
    Free(saar_ref);
    saar_ref = NULL;
    Free(delta_sa_ref);
    delta_sa_ref = NULL;
}

// PART 2: macros for wrappers
//
// Any reasonable C programmer should understand how these macros work,
// in light of Part 2. The number in the macro name is the number of 
// arguments sent to the C library function named 'cname'. See part 3
// for the handling of C library functions that return void.
//
// Note the check for NaN values within these macros; this is necessary
// because the GSW code does not check for NaN. (To be clear, in the
// GSW code, there are numerous comments about checking for NaN, but they
// refer to checking for GSW_INVALID_VALUE.)
#define W1(wname, cname, arg1, n, rval) \
void (wname)(double *(arg1), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        if (isnan((arg1)[i])) {\
            (rval)[i] = NA_REAL;\
        } else {\
            (rval)[i] = (cname)((arg1)[i]);\
            if ((rval)[i] == GSW_INVALID_VALUE) {\
                (rval)[i] = NA_REAL;\
            }\
        }\
    }\
}

#define W2(wname, cname, arg1, arg2, n, rval) \
void (wname)(double *(arg1), double *(arg2), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i])) {\
            (rval)[i] = NA_REAL;\
        } else {\
            (rval)[i] = (cname)((arg1)[i], (arg2)[i]);\
            if ((rval)[i] == GSW_INVALID_VALUE) {\
                (rval)[i] = NA_REAL;\
            }\
        }\
    }\
}

#define W3(wname, cname, arg1, arg2, arg3, n, rval) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i]) || isnan((arg3)[i])) {\
            (rval)[i] = NA_REAL;\
        } else {\
            (rval)[i] = (cname)((arg1)[i], (arg2)[i], (arg3)[i]);\
            if ((rval)[i] == GSW_INVALID_VALUE) {\
                (rval)[i] = NA_REAL;\
            }\
        }\
    }\
}

#define W4(wname, cname, arg1, arg2, arg3, arg4, n, rval) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), double *(arg4), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        /*Rprintf("args: %f %f %f %f %f; i=%d; n=%d\n", (arg1)[i], (arg2)[i], (arg3)[i], (arg4)[i], i, *(n));*/\
        /*Rprintf("invalid: %f\n", GSW_INVALID_VALUE);*/\
        /*if (isnan((arg1)[i])) Rprintf("arg1 is NaN\n");*/\
        /*if (isnan((arg2)[i])) Rprintf("arg2 is NaN\n");*/\
        /*if (isnan((arg3)[i])) Rprintf("arg3 is NaN\n");*/\
        /*if (isnan((arg4)[i])) Rprintf("arg4 is NaN\n");*/\
        if (isnan((arg1)[i]) || isnan((arg2)[i]) || isnan((arg3)[i]) || isnan((arg4)[i])) {\
            (rval)[i] = NA_REAL;\
        } else {\
            (rval)[i] = (cname)((arg1)[i], (arg2)[i], (arg3)[i], (arg4)[i]);\
            if ((rval)[i] == GSW_INVALID_VALUE) {\
                (rval)[i] = NA_REAL;\
            }\
        }\
    }\
}

#define W5(wname, cname, arg1, arg2, arg3, arg4, arg5, n, rval) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), double *(arg4), double *(arg5), int *(n), double *(rval))\
{\
    for (int i=0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i]) || isnan((arg3)[i]) || isnan((arg4)[i]) || isnn((arg5)[i])) {\
            (rval)[i] = NA_REAL;\
        } else {\
            (rval)[i] = (cname)((arg1)[i], (arg2)[i], (arg3)[i], (arg4)[i], (arg5)[i]);\
            if ((rval)[i] == GSW_INVALID_VALUE) {\
                (rval)[i] = NA_REAL;\
            }\
        }\
    }\
}


// PART 3: wrappers for functions that return a value. Wrapping is necessary 
// because the R function .C() cannot handle return values.
// See Part 3 for functios returning void.
//
// The case in the wrapper function names (first argument to macros)
// should match that in the TEOS-10 matlab functions, but note that case in
// the second argument should be lower-case, since this is used in the C
// library.
W3(wrap_gsw_adiabatic_lapse_rate_from_CT, gsw_adiabatic_lapse_rate_from_ct, SA, CT, p, n, rval)
W2(wrap_gsw_adiabatic_lapse_rate_ice, gsw_adiabatic_lapse_rate_ice, t, p, n, rval)
W3(wrap_gsw_alpha, gsw_alpha, SA, CT, p, n, rval)
W3(wrap_gsw_alpha_on_beta, gsw_alpha_on_beta, SA, CT, p, n, rval)
W3(wrap_gsw_alpha_wrt_t_exact, gsw_alpha_wrt_t_exact, SA, t, p, n, rval)
W2(wrap_gsw_alpha_wrt_t_ice, gsw_alpha_wrt_t_ice, t, p, n, rval)
W3(wrap_gsw_beta, gsw_beta, SA, CT, p, n, rval)
W3(wrap_gsw_beta_const_t_exact, gsw_beta_const_t_exact, SA, t, p, n, rval)
W3(wrap_gsw_cabbeling, gsw_cabbeling, SA, CT, p, n, rval)
W3(wrap_gsw_C_from_SP, gsw_c_from_sp, SP, t, p, n, rval)
W2(wrap_gsw_cp_ice, gsw_cp_ice, t, p, n, rval)
W3(wrap_gsw_cp_t_exact, gsw_cp_t_exact, SA, t, p, n, rval)
W3(wrap_gsw_CT_freezing, gsw_ct_freezing_exact, SA, p, saturation_fraction, n, rval)
W2(wrap_gsw_CT_from_pt, gsw_ct_from_pt, SA, pt, n, rval)
W3(wrap_gsw_CT_from_t, gsw_ct_from_t, SA, t, p, n, rval)
W4(wrap_gsw_deltaSA_from_SP, gsw_deltasa_from_sp, SP, p, longitude, latitude, n, rval)
W3(wrap_gsw_dynamic_enthalpy, gsw_dynamic_enthalpy, SA, CT, p, n, rval)
W3(wrap_gsw_enthalpy, gsw_enthalpy, SA, CT, p, n, rval)
W2(wrap_gsw_enthalpy_ice, gsw_enthalpy_ice, t, p, n, rval)
W3(wrap_gsw_enthalpy_t_exact, gsw_enthalpy_t_exact, SA, t, p, n, rval)
// gsw_entropy_from_CT is not in the C library.
W2(wrap_gsw_entropy_ice, gsw_entropy_ice, t, p, n, rval)
W3(wrap_gsw_entropy_from_t, gsw_entropy_from_t, SA, t, p, n, rval)
W2(wrap_gsw_grav, gsw_grav, latitude, p, n, rval)
W3(wrap_gsw_internal_energy, gsw_internal_energy, SA, CT, p, n, rval)
W2(wrap_gsw_internal_energy_ice, gsw_internal_energy_ice, t, p, n, rval)
W3(wrap_gsw_kappa, gsw_kappa, SA, CT, p, n, rval)
W2(wrap_gsw_kappa_const_t_ice, gsw_kappa_const_t_ice, t, p, n, rval)
W2(wrap_gsw_kappa_ice, gsw_kappa_ice, t, p, n, rval)
W3(wrap_gsw_kappa_t_exact, gsw_kappa_t_exact, SA, t, p, n, rval)
W2(wrap_gsw_latentheat_evap_CT, gsw_latentheat_evap_ct, SA, CT, n, rval)
W2(wrap_gsw_latentheat_evap_t, gsw_latentheat_evap_t, SA, t, n, rval)
W2(wrap_gsw_latentheat_melting, gsw_latentheat_melting, SA, p, n, rval)
W2(wrap_gsw_melting_ice_equilibrium_SA_CT_ratio, gsw_melting_ice_equilibrium_sa_ct_ratio, SA, p, n, rval)
W4(wrap_gsw_melting_ice_SA_CT_ratio, gsw_melting_ice_sa_ct_ratio, SA, CT, p, t_Ih, n, rval)
W2(wrap_gsw_pressure_coefficient_ice, gsw_pressure_coefficient_ice, t, p, n, rval)

// declare since it's not in the TEOS-10 C library yet, and was coded separately.
extern double gsw_p_from_z(double z, double latitude, double geo_strf_dyn_height, double sea_surface_geopotential);
W4(wrap_gsw_p_from_z, gsw_p_from_z, z, latitude, geo_strf_dyn_height, sea_surface_geopotential, n, rval)

W4(wrap_gsw_pot_rho_t_exact, gsw_pot_rho_t_exact, SA, t, p, p_ref, n, rval)
W3(wrap_gsw_pt0_from_t, gsw_pt0_from_t, SA, t, p, n, rval)
W2(wrap_gsw_pt_from_CT, gsw_pt_from_ct, SA, CT, n, rval)
W4(wrap_gsw_pt_from_t, gsw_pt_from_t, SA, t, p, p_ref, n, rval)
W3(wrap_gsw_rho, gsw_rho, SA, CT, p, n, rval)
W2(wrap_gsw_rho_ice, gsw_rho_ice, t, p, n, rval)
W3(wrap_gsw_rho_t_exact, gsw_rho_t_exact, SA, t, p, n, rval)
W3(wrap_gsw_SA_from_rho, gsw_sa_from_rho, rho, CT, p, n, rval)
W4(wrap_gsw_SA_from_SP, gsw_sa_from_sp, CT, p, longitude, latitude, n, rval)
W3(wrap_gsw_SA_from_SP_baltic, gsw_sa_from_sp_baltic, SP, longitude, latitude, n, rval)
W4(wrap_gsw_SA_from_Sstar, gsw_sa_from_sstar, Sstar, p, longitude, latitude, n, rval)
W1(wrap_gsw_SR_from_SP, gsw_sr_from_sp, SP, n, rval)
W2(wrap_gsw_sigma0, gsw_sigma0, SA, CT, n, rval)
W2(wrap_gsw_sigma1, gsw_sigma1, SA, CT, n, rval)
W2(wrap_gsw_sigma2, gsw_sigma2, SA, CT, n, rval)
W2(wrap_gsw_sigma3, gsw_sigma3, SA, CT, n, rval)
W2(wrap_gsw_sigma4, gsw_sigma4, SA, CT, n, rval)
W3(wrap_gsw_sound_speed, gsw_sound_speed, SA, t, p, n, rval)
W2(wrap_gsw_sound_speed_ice, gsw_sound_speed_ice, t, p, n, rval)
W3(wrap_gsw_sound_speed_t_exact, gsw_sound_speed_t_exact, SA, t, p, n, rval)
// gsw_specvol coded in R
W3(wrap_gsw_specvol_anom_standard, gsw_specvol_anom_standard, SA, CT, p, n, rval)
W2(wrap_gsw_specvol_ice, gsw_specvol_ice, t, p, n, rval)
W3(wrap_gsw_specvol_t_exact, gsw_specvol_t_exact, SA, t, p, n, rval)
W2(wrap_gsw_spiciness0, gsw_spiciness0, SA, CT, n, rval)
W2(wrap_gsw_spiciness1, gsw_spiciness1, SA, CT, n, rval)
W2(wrap_gsw_spiciness2, gsw_spiciness2, SA, CT, n, rval)
W3(wrap_gsw_SP_from_C, gsw_sp_from_c, C, t, p, n, rval)
W4(wrap_gsw_SP_from_SA, gsw_sp_from_sa, SA, p, longitude, latitude, n, rval)
W1(wrap_gsw_SP_from_SK, gsw_sp_from_sk, SK, n, rval)
W1(wrap_gsw_SP_from_SR, gsw_sp_from_sr, SR, n, rval)
W4(wrap_gsw_SP_from_Sstar, gsw_sp_from_sstar, Sstar, p, longitude, latitude, n, rval)
W4(wrap_gsw_Sstar_from_SA, gsw_sstar_from_sa, SA, p, longitude, latitude, n, rval)
W4(wrap_gsw_Sstar_from_SP, gsw_sstar_from_sp, SP, p, longitude, latitude, n, rval)
W3(wrap_gsw_t_freezing, gsw_t_freezing_exact, SA, p, saturation_fraction, n, rval)
//W3(wrap_gsw_t_freezing, gsw_t_freezing, SA, p, saturation_fraction, n, rval)
W3(wrap_gsw_t_from_CT, gsw_t_from_ct, SA, CT, p, n, rval)
W3(wrap_gsw_thermobaric, gsw_thermobaric, SA, CT, p, n, rval)
W2(wrap_gsw_z_from_p, gsw_z_from_p, p, lat, n, rval)

// PART 4
// Handle cases with multiple return values.
void wrap_gsw_CT_from_rho(double *rho, double *SA, double *p, int *n,
        double *CT, double *CT_multiple)
{
    for (int i=0; i < *(n); i++)
        gsw_ct_from_rho(rho[i], SA[i], p[i], &CT[i], &CT_multiple[i]);
}

void wrap_gsw_frazil_properties(double *SA_bulk, double *h_bulk, double *p, int *n,
        double *SA_final, double *CT_final, double *w_Ih_final)
{
    for (int i=0; i < *(n); i++)
        gsw_frazil_properties(SA_bulk[i], h_bulk[i], p[i],
                &SA_final[i], &CT_final[i], &w_Ih_final[i]);
}

void wrap_gsw_gibbs(int *ns, int *nt, int *np, double *SA, double *t, double *p, int *n, double *res)
{
    for (int i=0; i < *(n); i++)
        res[i] = gsw_gibbs(*ns, *nt, *np, SA[i], t[i], p[i]);
}

void wrap_gsw_gibbs_ice(int *nt, int *np, double *t, double *p, int *n, double *res)
{
    for (int i=0; i < *(n); i++)
        res[i] = gsw_gibbs_ice(*nt, *np, t[i], p[i]);
}

void wrap_gsw_ice_fraction_to_freeze_seawater(double *SA, double *CT, double *p, double *t_Ih, int *n,
        double *SA_freeze, double *CT_freeze, double *w_Ih)
{
    for (int i=0; i < *(n); i++)
        gsw_ice_fraction_to_freeze_seawater(SA[i], CT[i], p[i], t_Ih[i], &SA_freeze[i], &CT_freeze[i], &w_Ih[i]);
}

void wrap_gsw_IPV_vs_fNsquared_ratio(double *SA, double *CT, double *p, double *p_ref, int *n,
        double *IPV_vs_fNsquared_ratio, double *p_mid)
{
    gsw_ipv_vs_fnsquared_ratio(SA, CT, p, *p_ref, *n, IPV_vs_fNsquared_ratio, p_mid);
}

void wrap_gsw_rho_first_derivatives(double *SA, double *CT, double *p, int *n,
        double *drho_dsa, double *drho_dct, double *drho_dp)
{
    for (int i=0; i < *(n); i++)
        gsw_rho_first_derivatives(SA[i], CT[i], p[i], &drho_dsa[i], &drho_dct[i], &drho_dp[i]);
}

void wrap_gsw_rho_first_derivatives_wrt_enthalpy(double *SA, double *CT, double *p, int *n,
        double *rho_sa_wrt_h, double *rho_h)
{
    for (int i=0; i < *(n); i++)
        gsw_rho_first_derivatives_wrt_enthalpy(SA[i], CT[i], p[i], &rho_sa_wrt_h[i], &rho_h[i]);
}


// [SA_final, CT_final, w_Ih_final] = gsw_melting_ice_into_seawater(SA,CT,p,w_Ih,t_Ih)
void wrap_gsw_melting_ice_into_seawater(double *SA, double *CT, double *p, double *w_Ih, double *t_Ih,
        int *n, double *SA_final, double *CT_final, double *w_Ih_final)
{
    // see gswteos-10.h line 147 for declaration of gsw_melting_ice_into_seawater
    for (int i=0; i < *(n); i++) {
        gsw_melting_ice_into_seawater(SA[i], CT[i], p[i], w_Ih[i], t_Ih[i],
                &SA_final[i], &CT_final[i], &w_Ih_final[i]);
    }
}

void wrap_gsw_Nsquared(double *SA, double *CT, double *p, double *latitude, int *n, double *n2, double *p_mid)
{
    extern void gsw_nsquared(double *sa, double *ct, double *p, double *latitude, int nz, double *n2, double *p_mid);
    gsw_nsquared(SA, CT, p, latitude, *n, n2, p_mid);
}

void wrap_gsw_rho_alpha_beta(double *SA, double *CT, double *p, int *n,
        double *rho, double *alpha, double *beta)
{
    for (int i=0; i < *(n); i++)
        gsw_rho_alpha_beta(SA[i], CT[i], p[i], &rho[i], &alpha[i], &beta[i]);
}

void wrap_gsw_SAAR(double *p, double *longitude, double *latitude, int *n, double *saar, int *inocean)
{
    for (int i=0; i < *(n); i++) {
        saar[i] = gsw_saar(p[i], longitude[i], latitude[i]);
        inocean[i] = saar[i]==0?0:1;
    }
}

void wrap_gsw_Turner_Rsubrho(double *SA, double *CT, double *p, int *n, double *Tu, double *Rsubrho, double *p_mid)
{
    extern void gsw_turner_rsubrho(double *sa, double *ct, double *p, int nz, double *Tu, double *Rsubrho, double *p_mid);
    gsw_turner_rsubrho(SA, CT, p, *n, Tu, Rsubrho, p_mid);
}

