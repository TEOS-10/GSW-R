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

// Wrappers are used because the R function .C() cannot handle
// return values. A few notes may help the reader.
//
// 1. Macros (e.g. W11, etc) are used for most cases,
// but some, such as gsw_gibbs(), have argument lists that are not
// like the majority of the other functions, so no macros are used in that
// case.
//
// 2. Macros are only used for C functions that return either a single
// scalar or a set of scalars. They are not used for C functions that fill
// up vectors (e.g. the N^2 calculation).
//
// 3. The case in the wrapper function names (first argument to macros)
// should match that in the TEOS-10 matlab functions, but that case in
// the second argument should be lower-case, since this is the convention
// of the underlying C library that does all the actual calculations.
//
// 4. We convert NaN values to the GSW convention for missing values,
// which is what is checked for whenever C code checks for 
// GSW_INVALID_VALUE.
//
// 5. Macro names reflect the number of input and output values. The first
// digit after "W" is the number of input values (not counting "n", which
// is the vector length of each of these input values) and the second digit
// is the number of output values. Thus, for example, W32 takes three input
// values (plus "n") and returns 2 output values.

#define W11(wname, cname, arg1, n, rval) \
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

// 2 input values, 1 output value
#define W21(wname, cname, arg1, arg2, n, rval) \
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

// 2 input values, 2 output values
#define W22(wname, cname, arg1, arg2, n, rval1, rval2) \
void (wname)(double *(arg1), double *(arg2), int *(n), double *(rval1), double *(rval2))\
{\
    for (int i = 0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i])) {\
            (rval1)[i] = NA_REAL;\
            (rval2)[i] = NA_REAL;\
        } else {\
            (cname)((arg1)[i], (arg2)[i], &(rval1)[i], &(rval2)[i]);\
            if ((rval1)[i] == GSW_INVALID_VALUE) (rval1)[i] = NA_REAL;\
            if ((rval2)[i] == GSW_INVALID_VALUE) (rval2)[i] = NA_REAL;\
        }\
    }\
}

// 2 input values, 3 output values
#define W23(wname, cname, arg1, arg2, n, rval1, rval2, rval3) \
void (wname)(double *(arg1), double *(arg2), int *(n), double *(rval1), double *(rval2), double *(rval3))\
{\
    for (int i = 0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i])) {\
            (rval1)[i] = NA_REAL;\
            (rval2)[i] = NA_REAL;\
            (rval3)[i] = NA_REAL;\
        } else {\
            (cname)((arg1)[i], (arg2)[i], &(rval1)[i], &(rval2)[i], &(rval3)[i]);\
            if ((rval1)[i] == GSW_INVALID_VALUE) (rval1)[i] = NA_REAL;\
            if ((rval2)[i] == GSW_INVALID_VALUE) (rval2)[i] = NA_REAL;\
            if ((rval3)[i] == GSW_INVALID_VALUE) (rval3)[i] = NA_REAL;\
        }\
    }\
}


// 3 input values, 1 output value
#define W31(wname, cname, arg1, arg2, arg3, n, rval) \
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

// 3 input values, 2 output values
#define W32(wname, cname, arg1, arg2, arg3, n, rval1, rval2) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), int *(n), double *(rval1), double *(rval2))\
{\
    for (int i = 0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i]) || isnan((arg3)[i])) {\
            (rval1)[i] = NA_REAL;\
            (rval2)[i] = NA_REAL;\
        } else {\
            (cname)((arg1)[i], (arg2)[i], (arg3)[i], &(rval1)[i], &(rval2)[i]);\
            if ((rval1)[i] == GSW_INVALID_VALUE) (rval1)[i] = NA_REAL;\
            if ((rval2)[i] == GSW_INVALID_VALUE) (rval2)[i] = NA_REAL;\
        }\
    }\
}

// 3 input values, 3 output values
#define W33(wname, cname, arg1, arg2, arg3, n, rval1, rval2, rval3) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), int *(n), double *(rval1), double *(rval2), double *(rval3))\
{\
    for (int i = 0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i]) || isnan((arg3)[i])) {\
            (rval1)[i] = NA_REAL;\
            (rval2)[i] = NA_REAL;\
            (rval3)[i] = NA_REAL;\
        } else {\
            (cname)((arg1)[i], (arg2)[i], (arg3)[i], &(rval1)[i], &(rval2)[i], &(rval3)[i]);\
            if ((rval1)[i] == GSW_INVALID_VALUE) (rval1)[i] = NA_REAL;\
            if ((rval2)[i] == GSW_INVALID_VALUE) (rval2)[i] = NA_REAL;\
            if ((rval3)[i] == GSW_INVALID_VALUE) (rval3)[i] = NA_REAL;\
        }\
    }\
}

// 3 input values, 5 output values
#define W35(wname, cname, arg1, arg2, arg3, n, rval1, rval2, rval3, rval4, rval5) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), int *(n), double *(rval1), double *(rval2), double *(rval3), double *(rval4), double *(rval5))\
{\
    for (int i = 0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i]) || isnan((arg3)[i])) {\
            (rval1)[i] = NA_REAL;\
            (rval2)[i] = NA_REAL;\
            (rval3)[i] = NA_REAL;\
            (rval4)[i] = NA_REAL;\
            (rval5)[i] = NA_REAL;\
        } else {\
            (cname)((arg1)[i], (arg2)[i], (arg3)[i], &(rval1)[i], &(rval2)[i], &(rval3)[i], &(rval4)[i], &(rval5)[i]);\
            if ((rval1)[i] == GSW_INVALID_VALUE) (rval1)[i] = NA_REAL;\
            if ((rval2)[i] == GSW_INVALID_VALUE) (rval2)[i] = NA_REAL;\
            if ((rval3)[i] == GSW_INVALID_VALUE) (rval3)[i] = NA_REAL;\
            if ((rval4)[i] == GSW_INVALID_VALUE) (rval4)[i] = NA_REAL;\
            if ((rval5)[i] == GSW_INVALID_VALUE) (rval5)[i] = NA_REAL;\
        }\
    }\
}


// 4 input parameters, 1 output value
#define W41(wname, cname, arg1, arg2, arg3, arg4, n, rval) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), double *(arg4), int *(n), double *(rval))\
{\
    for (int i = 0; i < *(n); i++) {\
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

// 4 input parameters, 3 output values
#define W43(wname, cname, arg1, arg2, arg3, arg4, n, rval1, rval2, rval3) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), double *(arg4), int *(n), double *(rval1), double *(rval2), double *(rval3))\
{\
    for (int i = 0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i]) || isnan((arg3)[i]) || isnan((arg4)[i])) {\
            (rval1)[i] = NA_REAL;\
            (rval2)[i] = NA_REAL;\
            (rval3)[i] = NA_REAL;\
        } else {\
            (cname)((arg1)[i], (arg2)[i], (arg3)[i], (arg4)[i], &(rval1)[i], &(rval2)[i], &(rval3)[i]);\
            if ((rval1)[i] == GSW_INVALID_VALUE) (rval1)[i] = NA_REAL;\
            if ((rval2)[i] == GSW_INVALID_VALUE) (rval2)[i] = NA_REAL;\
            if ((rval3)[i] == GSW_INVALID_VALUE) (rval3)[i] = NA_REAL;\
        }\
    }\
}

// 5 input parameters, 3 output values
#define W53(wname, cname, arg1, arg2, arg3, arg4, arg5, n, rval1, rval2, rval3) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), double *(arg4), double *(arg5), int *(n), double *(rval1), double *(rval2), double *(rval3))\
{\
    for (int i = 0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i]) || isnan((arg3)[i]) || isnan((arg4)[i]) || isnan((arg5)[i])) {\
            (rval1)[i] = NA_REAL;\
            (rval2)[i] = NA_REAL;\
            (rval3)[i] = NA_REAL;\
        } else {\
            (cname)((arg1)[i], (arg2)[i], (arg3)[i], (arg4)[i], (arg5)[i], &(rval1)[i], &(rval2)[i], &(rval3)[i]);\
            if ((rval1)[i] == GSW_INVALID_VALUE) (rval1)[i] = NA_REAL;\
            if ((rval2)[i] == GSW_INVALID_VALUE) (rval2)[i] = NA_REAL;\
            if ((rval3)[i] == GSW_INVALID_VALUE) (rval3)[i] = NA_REAL;\
        }\
    }\
}

// 6 input parameters, 2 output values
#define W63(wname, cname, arg1, arg2, arg3, arg4, arg5, arg6, n, rval1, rval2) \
void (wname)(double *(arg1), double *(arg2), double *(arg3), double *(arg4), double *(arg5), double *(arg6), int *(n), double *(rval1), double *(rval2))\
{\
    for (int i = 0; i < *(n); i++) {\
        if (isnan((arg1)[i]) || isnan((arg2)[i]) || isnan((arg3)[i]) || isnan((arg4)[i]) || isnan((arg5)[i]) || isnan((arg6)[i])) {\
            (rval1)[i] = NA_REAL;\
            (rval2)[i] = NA_REAL;\
        } else {\
            (cname)((arg1)[i], (arg2)[i], (arg3)[i], (arg4)[i], (arg5)[i], (arg6)[i], &(rval1)[i], &(rval2)[i]);\
            if ((rval1)[i] == GSW_INVALID_VALUE) (rval1)[i] = NA_REAL;\
            if ((rval2)[i] == GSW_INVALID_VALUE) (rval2)[i] = NA_REAL;\
        }\
    }\
}
 
 
W31(wrap_gsw_adiabatic_lapse_rate_from_CT, gsw_adiabatic_lapse_rate_from_ct, SA, CT, p, n, rval)
W21(wrap_gsw_adiabatic_lapse_rate_ice, gsw_adiabatic_lapse_rate_ice, t, p, n, rval)
W31(wrap_gsw_alpha, gsw_alpha, SA, CT, p, n, rval)
W31(wrap_gsw_alpha_on_beta, gsw_alpha_on_beta, SA, CT, p, n, rval)
W31(wrap_gsw_alpha_wrt_t_exact, gsw_alpha_wrt_t_exact, SA, t, p, n, rval)
W21(wrap_gsw_alpha_wrt_t_ice, gsw_alpha_wrt_t_ice, t, p, n, rval)
W31(wrap_gsw_beta, gsw_beta, SA, CT, p, n, rval)
W31(wrap_gsw_beta_const_t_exact, gsw_beta_const_t_exact, SA, t, p, n, rval)
W31(wrap_gsw_cabbeling, gsw_cabbeling, SA, CT, p, n, rval)
W31(wrap_gsw_C_from_SP, gsw_c_from_sp, SP, t, p, n, rval)
W21(wrap_gsw_chem_potential_water_ice, gsw_chem_potential_water_ice, t, p, n, rval)
W31(wrap_gsw_chem_potential_water_t_exact, gsw_chem_potential_water_t_exact, SA, t, p, n, rval)
W21(wrap_gsw_cp_ice, gsw_cp_ice, t, p, n, rval)
W31(wrap_gsw_cp_t_exact, gsw_cp_t_exact, SA, t, p, n, rval)
W22(wrap_gsw_CT_first_derivatives, gsw_ct_first_derivatives, SA, pt, n, CT_SA, CT_pt)
W33(wrap_gsw_CT_first_derivatives_wrt_t_exact, gsw_ct_first_derivatives_wrt_t_exact, SA, t, p, n, CT_SA_wrt_t, CT_t_wrt_t, CT_p_wrt_t)
W31(wrap_gsw_CT_freezing, gsw_ct_freezing, SA, p, saturation_fraction, n, rval)
//W31(wrap_gsw_CT_freezing_exact, gsw_ct_freezing_exact, SA, p, saturation_fraction, n, rval)
W32(wrap_gsw_CT_freezing_first_derivatives, gsw_ct_freezing_first_derivatives, SA, p, saturation_fraction, n, CTfreezing_SA, CTfreezing_p)
W32(wrap_gsw_CT_freezing_first_derivatives_poly,gsw_ct_freezing_first_derivatives_poly,SA,p,saturation_fraction,n,CTfreezing_SA,CTfreezing_p)
W31(wrap_gsw_CT_from_enthalpy, gsw_ct_from_enthalpy, SA, h, p, n, rval)
W21(wrap_gsw_CT_from_entropy, gsw_ct_from_entropy, SA, entropy, n, rval)
W21(wrap_gsw_CT_from_pt, gsw_ct_from_pt, SA, pt, n, rval)
W32(wrap_gsw_CT_from_rho, gsw_ct_from_rho, rho, SA, p, n, CT, CT_multiple)
W31(wrap_gsw_CT_from_t, gsw_ct_from_t, SA, t, p, n, rval)
W21(wrap_gsw_CT_maxdensity, gsw_ct_maxdensity, SA, p, n, rval)
W23(wrap_gsw_CT_second_derivatives, gsw_ct_second_derivatives, SA, pt, n, CT_SA_SA, CT_SA_pt, CT_pt_pt)
W41(wrap_gsw_deltaSA_from_SP, gsw_deltasa_from_sp, SP, p, longitude, latitude, n, rval)
W31(wrap_gsw_dilution_coefficient_t_exact, gsw_dilution_coefficient_t_exact, SA, t, p, n, rval)
W31(wrap_gsw_dynamic_enthalpy, gsw_dynamic_enthalpy, SA, CT, p, n, rval)
W31(wrap_gsw_enthalpy, gsw_enthalpy, SA, CT, p, n, rval)
W31(wrap_gsw_enthalpy_ct_exact, gsw_enthalpy_ct_exact, SA, t, p, n, rval)
W41(wrap_gsw_enthalpy_diff, gsw_enthalpy_diff, SA, CT, p_shallow, p_deep, n, rval)
W32(wrap_gsw_enthalpy_first_derivatives, gsw_enthalpy_first_derivatives, SA, CT, p, n, h_SA, h_CT)
W32(wrap_gsw_enthalpy_first_derivatives_CT_exact, gsw_enthalpy_first_derivatives_ct_exact, SA, CT, p, n, h_SA, h_CT)
W21(wrap_gsw_enthalpy_ice, gsw_enthalpy_ice, t, p, n, rval)
W33(wrap_gsw_enthalpy_second_derivatives, gsw_enthalpy_second_derivatives, SA, CT, p, n, h_SA_SA, h_SA_CT, h_CT_CT)
W33(wrap_gsw_enthalpy_second_derivatives_CT_exact, gsw_enthalpy_second_derivatives_ct_exact, SA, CT, p, n, h_SA_SA, h_SA_CT, h_CT_CT)
W31(wrap_gsw_enthalpy_t_exact, gsw_enthalpy_t_exact, SA, t, p, n, rval)
W22(wrap_gsw_entropy_first_derivatives, gsw_entropy_first_derivatives, SA, CT, n, eta_SA, eta_CT)
// gsw_entropy_from_CT is not in the C library.
W21(wrap_gsw_entropy_ice, gsw_entropy_ice, t, p, n, rval)
W21(wrap_gsw_entropy_from_pt, gsw_entropy_from_pt, SA, pt, n, rval)
W31(wrap_gsw_entropy_from_t, gsw_entropy_from_t, SA, t, p, n, rval)
W23(wrap_gsw_entropy_second_derivatives, gsw_entropy_second_derivatives, SA, CT, n, h_SA_SA, h_SA_CT, h_CT_CT)
W31(wrap_gsw_Fdelta, gsw_fdelta, p, longitude, latitude, n, rval)
W33(wrap_gsw_frazil_properties, gsw_frazil_properties, SA_bulk, h_bulk, p, n, SA_final, CT_final, w_Ih_final)
W33(wrap_gsw_frazil_properties_potential, gsw_frazil_properties_potential, SA_bulk, h_pot_bulk, p, n, SA_final, CT_final, w_Ih_final)
W33(wrap_gsw_frazil_properties_potential_poly, gsw_frazil_properties_potential_poly, SA_bulk, h_pot_bulk, p, n, SA_final, CT_final, w_Ih_final)
W33(wrap_gsw_frazil_ratios_adiabatic, gsw_frazil_ratios_adiabatic, SA, p, w_Ih, n, dSA_dCT_frazil, dSA_dP_frazil, dCT_dP_frazil)
W33(wrap_gsw_frazil_ratios_adiabatic_poly, gsw_frazil_ratios_adiabatic_poly, SA, p, w_Ih, n, dSA_dCT_frazil, dSA_dP_frazil, dCT_dP_frazil)

void wrap_gsw_geo_strf_dyn_height(double *SA, double *CT, double *p, double *p_ref, int *n, double *dyn_height)
{
    gsw_geo_strf_dyn_height(SA, CT, p, *p_ref, *n, dyn_height);
}

void wrap_gsw_geo_strf_dyn_height_pc(double *SA, double *CT, double *delta_p, int *n, double *dyn_height, double *p_mid)
{
    gsw_geo_strf_dyn_height_pc(SA, CT, delta_p, *n, dyn_height, p_mid);
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
W21(wrap_gsw_grav, gsw_grav, latitude, p, n, rval)
W21(wrap_gsw_Helmholtz_energy_ice, gsw_helmholtz_energy_ice, t, p, n, rval)
W11(wrap_gsw_hill_ratio_at_sp2, gsw_hill_ratio_at_sp2, t, n, rval)
W43(wrap_gsw_ice_fraction_to_freeze_seawater, gsw_ice_fraction_to_freeze_seawater, SA, CT, p, t_Ih, n, SA_freeze, CT_freeze, w_Ih)
W31(wrap_gsw_internal_energy, gsw_internal_energy, SA, CT, p, n, rval)
W21(wrap_gsw_internal_energy_ice, gsw_internal_energy_ice, t, p, n, rval)
void wrap_gsw_IPV_vs_fNsquared_ratio(double *SA, double *CT, double *p, double *p_ref, int *n,
        double *IPV_vs_fNsquared_ratio, double *p_mid)
{
    gsw_ipv_vs_fnsquared_ratio(SA, CT, p, *p_ref, *n, IPV_vs_fNsquared_ratio, p_mid);
}
W31(wrap_gsw_kappa, gsw_kappa, SA, CT, p, n, rval)
W21(wrap_gsw_kappa_const_t_ice, gsw_kappa_const_t_ice, t, p, n, rval)
W21(wrap_gsw_kappa_ice, gsw_kappa_ice, t, p, n, rval)
W31(wrap_gsw_kappa_t_exact, gsw_kappa_t_exact, SA, t, p, n, rval)
W21(wrap_gsw_latentheat_evap_CT, gsw_latentheat_evap_ct, SA, CT, n, rval)
W21(wrap_gsw_latentheat_evap_t, gsw_latentheat_evap_t, SA, t, n, rval)
W21(wrap_gsw_latentheat_melting, gsw_latentheat_melting, SA, p, n, rval)
W21(wrap_gsw_melting_ice_equilibrium_SA_CT_ratio, gsw_melting_ice_equilibrium_sa_ct_ratio, SA, p, n, rval)
W21(wrap_gsw_melting_ice_equilibrium_SA_CT_ratio_poly, gsw_melting_ice_equilibrium_sa_ct_ratio_poly, SA, p, n, rval)
W53(wrap_gsw_melting_ice_into_seawater, gsw_melting_ice_into_seawater, SA, CT, p, w_Ih, t_Ih, n, SA_final, CT_final, w_Ih_final)
W41(wrap_gsw_melting_ice_SA_CT_ratio, gsw_melting_ice_sa_ct_ratio, SA, CT, p, t_Ih, n, rval)
W41(wrap_gsw_melting_ice_SA_CT_ratio_poly, gsw_melting_ice_sa_ct_ratio_poly, SA, CT, p, t_Ih, n, rval)
W63(wrap_gsw_melting_seaice_into_seawater, gsw_melting_seaice_into_seawater, SA, CT, p, w_seaice, SA_seaice, t_seaice, n, SA_final, CT_final)
void wrap_gsw_Nsquared(double *SA, double *CT, double *p, double *latitude, int *n, double *n2, double *p_mid)
{
    extern void gsw_nsquared(double *sa, double *ct, double *p, double *latitude, int nz, double *n2, double *p_mid);
    gsw_nsquared(SA, CT, p, latitude, *n, n2, p_mid);
}
W11(wrap_gsw_pot_enthalpy_from_pt_ice, gsw_pot_enthalpy_from_pt_ice, pt_ice, n, rval)
W11(wrap_gsw_pot_enthalpy_from_pt_ice_poly, gsw_pot_enthalpy_from_pt_ice_poly, pt_ice, n, rval)
//W31(wrap_gsw_pot_enthalpy_ice_freezing, gsw_pot_enthalpy_ice_freezing, SA, p, saturaion_fraction, n, rval)
W21(wrap_gsw_pot_enthalpy_ice_freezing, gsw_pot_enthalpy_ice_freezing, SA, p, n, rval)
W21(wrap_gsw_pot_enthalpy_ice_freezing_poly, gsw_pot_enthalpy_ice_freezing_poly, SA, p, n, rval)
W22(wrap_gsw_pot_enthalpy_ice_freezing_first_derivatives, gsw_pot_enthalpy_ice_freezing_first_derivatives, SA, p, n, pot_enthalpy_ice_freezing_SA, pot_enthalpy_ice_freezing_p)
W22(wrap_gsw_pot_enthalpy_ice_freezing_first_derivatives_poly, gsw_pot_enthalpy_ice_freezing_first_derivatives_poly, SA, p, n, pot_enthalpy_ice_freezing_SA, pot_enthalpy_ice_freezing_p)
W21(wrap_gsw_pressure_coefficient_ice, gsw_pressure_coefficient_ice, t, p, n, rval)
W31(wrap_gsw_pressure_freezing_CT, gsw_pressure_freezing_ct, SA, CT, saturation_fraction, n, rval)
// The next line is necessary because gsw_p_from_z() is not in the TEOS-10 C library yet.
extern double gsw_p_from_z(double z, double latitude);
W21(wrap_gsw_p_from_z, gsw_p_from_z, z, latitude, n, rval)
W41(wrap_gsw_pot_rho_t_exact, gsw_pot_rho_t_exact, SA, t, p, p_ref, n, rval)
W31(wrap_gsw_pt0_from_t, gsw_pt0_from_t, SA, t, p, n, rval)
W21(wrap_gsw_pt0_from_t_ice, gsw_pt0_from_t_ice, t, p, n, rval)
W22(wrap_gsw_pt_first_derivatives, gsw_pt_first_derivatives, SA, CT, n, pt_SA, pt_CT)
W21(wrap_gsw_pt_from_CT, gsw_pt_from_ct, SA, CT, n, rval)
W21(wrap_gsw_pt_from_entropy, gsw_pt_from_entropy, SA, entropy, n, rval)
W11(wrap_gsw_pt_from_pot_enthalpy_ice, gsw_pt_from_pot_enthalpy_ice, pot_enthalpy_ice, n, rval)
W11(wrap_gsw_pt_from_pot_enthalpy_ice_poly, gsw_pt_from_pot_enthalpy_ice_poly, pot_enthalpy_ice, n, rval)
W41(wrap_gsw_pt_from_t, gsw_pt_from_t, SA, t, p, p_ref, n, rval)
W31(wrap_gsw_pt_from_t_ice, gsw_pt_from_t_ice, t, p, p_ref, n, rval)
W23(wrap_gsw_pt_second_derivatives, gsw_pt_second_derivatives, SA, CT, n, pt_SA_SA, pt_SA_CT, pt_CT_CT)
W31(wrap_gsw_rho, gsw_rho, SA, CT, p, n, rval)
W33(wrap_gsw_rho_alpha_beta, gsw_rho_alpha_beta, SA, CT, p, n, rho, alpha, beta)
W33(wrap_gsw_rho_first_derivatives, gsw_rho_first_derivatives, SA, CT, p, n, drho_dsa, drho_dct, drho_dp)
W32(wrap_gsw_rho_first_derivatives_wrt_enthalpy, gsw_rho_first_derivatives_wrt_enthalpy, SA, CT, p, n, rho_sa_wrt_h, rho_h)
W21(wrap_gsw_rho_ice, gsw_rho_ice, t, p, n, rval)
W35(wrap_gsw_rho_second_derivatives, gsw_rho_second_derivatives, SA, CT, p, n, rho_SA_SA, rho_SA_CT, rho_CT_CT, rho_SA_p, rho_CT_p)
W33(wrap_gsw_rho_second_derivatives_wrt_enthalpy, gsw_rho_second_derivatives_wrt_enthalpy, SA, CT, p, n, rho_SA_SA, rho_SA_h, rho_h_h)
W31(wrap_gsw_rho_t_exact, gsw_rho_t_exact, SA, t, p, n, rval)
void wrap_gsw_SAAR(double *p, double *longitude, double *latitude, int *n, double *saar, int *inocean)
{
    for (int i=0; i < *(n); i++) {
        saar[i] = gsw_saar(p[i], longitude[i], latitude[i]);
        inocean[i] = saar[i]==0?0:1;
    }
}
W31(wrap_gsw_SA_freezing_from_CT, gsw_sa_freezing_from_ct, CT, p, saturation_fraction, n, rval)
W31(wrap_gsw_SA_freezing_from_CT_poly, gsw_sa_freezing_from_ct_poly, CT, p, saturation_fraction, n, rval)
W31(wrap_gsw_SA_freezing_from_t, gsw_sa_freezing_from_t, t, p, saturation_fraction, n, rval)
W31(wrap_gsw_SA_freezing_from_t_poly, gsw_sa_freezing_from_t_poly, t, p, saturation_fraction, n, rval)
W31(wrap_gsw_SA_from_rho, gsw_sa_from_rho, rho, CT, p, n, rval)
W41(wrap_gsw_SA_from_SP, gsw_sa_from_sp, CT, p, longitude, latitude, n, rval)
W31(wrap_gsw_SA_from_SP_Baltic, gsw_sa_from_sp_baltic, SP, longitude, latitude, n, rval)
W41(wrap_gsw_SA_from_Sstar, gsw_sa_from_sstar, Sstar, p, longitude, latitude, n, rval)
W53(wrap_gsw_seaice_fraction_to_freeze_seawater, gsw_seaice_fraction_to_freeze_seawater, SA, CT, p, SA_seaice, t_seaice, n, SA_freeze, CT_freeze, w_seaice)
W11(wrap_gsw_SR_from_SP, gsw_sr_from_sp, SP, n, rval)
W21(wrap_gsw_sigma0, gsw_sigma0, SA, CT, n, rval)
W21(wrap_gsw_sigma1, gsw_sigma1, SA, CT, n, rval)
W21(wrap_gsw_sigma2, gsw_sigma2, SA, CT, n, rval)
W21(wrap_gsw_sigma3, gsw_sigma3, SA, CT, n, rval)
W21(wrap_gsw_sigma4, gsw_sigma4, SA, CT, n, rval)
W31(wrap_gsw_sound_speed, gsw_sound_speed, SA, t, p, n, rval)
W21(wrap_gsw_sound_speed_ice, gsw_sound_speed_ice, t, p, n, rval)
W31(wrap_gsw_sound_speed_t_exact, gsw_sound_speed_t_exact, SA, t, p, n, rval)
W33(wrap_gsw_specvol_alpha_beta, gsw_specvol_alpha_beta, SA, CT, p, n, specvol, alpha, beta)
// gsw_specvol coded in R (??)
W31(wrap_gsw_specvol_anom_standard, gsw_specvol_anom_standard, SA, CT, p, n, rval)
W33(wrap_gsw_specvol_first_derivatives, gsw_specvol_first_derivatives, SA, CT, p, n, v_SA, v_CT, v_p)
W32(wrap_gsw_specvol_first_derivatives_wrt_enthalpy, gsw_specvol_first_derivatives_wrt_enthalpy, SA, CT, p, n, v_SA, v_h)
W21(wrap_gsw_specvol_ice, gsw_specvol_ice, t, p, n, rval)

W35(wrap_gsw_specvol_second_derivatives, gsw_specvol_second_derivatives, SA, CT, p, n, specvol_SA_SA, specvol_SA_CT, specvol_CT_CT, specvol_SA_p, specvol_CT_p)
W33(wrap_gsw_specvol_second_derivatives_wrt_enthalpy, gsw_specvol_second_derivatives_wrt_enthalpy, SA, CT, p, n, specvol_SA_SA, specvol_SA_h, specvol_h_h)

W31(wrap_gsw_specvol_t_exact, gsw_specvol_t_exact, SA, t, p, n, rval)
W21(wrap_gsw_spiciness0, gsw_spiciness0, SA, CT, n, rval)
W21(wrap_gsw_spiciness1, gsw_spiciness1, SA, CT, n, rval)
W21(wrap_gsw_spiciness2, gsw_spiciness2, SA, CT, n, rval)
W31(wrap_gsw_SP_from_C, gsw_sp_from_c, C, t, p, n, rval)
W41(wrap_gsw_SP_from_SA, gsw_sp_from_sa, SA, p, longitude, latitude, n, rval)
W11(wrap_gsw_SP_from_SK, gsw_sp_from_sk, SK, n, rval)
W11(wrap_gsw_SP_from_SR, gsw_sp_from_sr, SR, n, rval)
W41(wrap_gsw_SP_from_Sstar, gsw_sp_from_sstar, Sstar, p, longitude, latitude, n, rval)
W41(wrap_gsw_Sstar_from_SA, gsw_sstar_from_sa, SA, p, longitude, latitude, n, rval)
W41(wrap_gsw_Sstar_from_SP, gsw_sstar_from_sp, SP, p, longitude, latitude, n, rval)
W31(wrap_gsw_t_deriv_chem_potential_water_t_exact, gsw_t_deriv_chem_potential_water_t_exact, SA, t, p, n, rval)
W31(wrap_gsw_t_freezing, gsw_t_freezing, SA, p, saturation_fraction, n, rval)
W32(wrap_gsw_t_freezing_first_derivatives, gsw_t_freezing_first_derivatives, SA, p, saturation_fraction, n, tfreezing_SA, tfreezing_p)
W32(wrap_gsw_t_freezing_first_derivatives_poly, gsw_t_freezing_first_derivatives_poly, SA, p, saturation_fraction, n, tfreezing_SA, tfreezing_p)
W31(wrap_gsw_t_from_CT, gsw_t_from_ct, SA, CT, p, n, rval)
W21(wrap_gsw_t_from_pt0_ice, gsw_t_from_pt0_ice, pt0_ice, p, n, rval)
W31(wrap_gsw_thermobaric, gsw_thermobaric, SA, CT, p, n, rval)
// gsw_turner_rsubrho() works on an *entire* profile; it should not be wrapped in W33.
void wrap_gsw_Turner_Rsubrho(double *SA, double *CT, double *p, int *n, double *Tu, double *Rsubrho, double *p_mid)
{
    extern void gsw_turner_rsubrho(double *sa, double *ct, double *p, int nz, double *Tu, double *Rsubrho, double *p_mid);
    gsw_turner_rsubrho(SA, CT, p, *n, Tu, Rsubrho, p_mid);
}
W21(wrap_gsw_z_from_p, gsw_z_from_p, p, lat, n, rval)

