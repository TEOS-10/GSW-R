#include "gswteos-10.h"
double
gsw_p_from_z(double z, double latitude,
        double geo_strf_dyn_height, double sea_surface_geopotental)
{
    // Comments are from the gsw_matlab_v3_04 file gsw_p_from_z.m, which is the 
    // source for this transliteration into C by Dan Kelley. Most likely a new
    // version of the C library will replace this at some point in 2015.
    double db2Pa = 1e4; 
    double gamma = 2.26e-7;
    // If the graviational acceleration were to be regarded as 
    // being depth-independent, which is often the case in 
    // ocean models, then gamma would be set to be zero here,
    // and the code below works perfectly well.  
    double deg2rad = M_PI / 180.0;
    double X = sin(latitude * deg2rad);
    double sin2 = X * X;
    double gs = 9.780327*(1.0 + (5.2792e-3 + (2.32e-5 * sin2)) * sin2);

    // get the first estimate of p from Saunders (1981)
    double c1 = 5.25e-3 * sin2 + 5.92e-3;
    double p = -2 * z / ((1 - c1) + sqrt((1 - c1) * (1 - c1) + 8.84e-6 * z));
    // end of the first estimate from Saunders (1981)

    double df_dp = db2Pa * gsw_specvol_sso_0_p(p);
    // initial value of the derivative of f

    double f = gsw_enthalpy_sso_0_p(p) + gs * (z - 0.5*gamma*(z * z)) 
        - (geo_strf_dyn_height + sea_surface_geopotental);
    double p_old = p;
    p = p_old - f / df_dp;
    double p_mid = 0.5 * (p + p_old);
    df_dp = db2Pa * gsw_specvol_sso_0_p(p_mid);
    p = p_old - f / df_dp;
    return(p);
}

