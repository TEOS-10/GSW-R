## Please see ../README_developer.md for the scheme used in adding functions
## here. Generally the functions will be added to Part 4.


## PART 1: document the package and the 'gsw' dataset

#' R implementation of the Thermodynamic Equation Of Seawater - 2010 (TEOS-10)
#'
#' @description
#' Provides an R interface to the TEOS-10 / GSW (Gibbs Sea Water) library,
#' partly for use by the \code{oce} package (see \url{http://dankelley.github.io/oce})
#' and partly for general use. It is assumed that users are familiar with
#' the science and methodology of GSW, and that the package vignette
#' (obtained by typing \code{vignette("gsw")} in an R window) provides 
#' enough orientation to get users started with the \code{gsw} functions.
#'
#' @details
#' \code{gsw} was developed using open-source methodologies, on
#' the GitHub site (\url{https://github.com/TEOS-10/GSW-R}), which
#' is part of a set of sites dedicated to GSW formulations in 
#' various languages.
#'
#' The \code{gsw} system is to link R functions with the C version of
#' the TEOS-10 library.  The R function names are chosen to match
#' those of the Matlab version of GSW, and the function arguments
#' also match with one exception: in \code{gsw}, longitude
#' and latitude are indicated with their full names, whereas in
#' Matlab they are indicated with \code{long} and \code{lat};
#' since R permits abbreviated function arguments, the shortened
#' names can be used in \code{gsw} as well.
#'
#' The documentation for the \code{gsw} functions focuses mainly
#' on the arguments and return values, relying on links to the
#' TEOS-10 webpages for details.
#' 
#' See \url{http://www.teos-10.org/pubs/gsw/html/gsw_contents.html}
#' for a list of the TEOS-10 functions and 
#' \url{http://teos-10.github.io/GSW-R/documentation.html} for a list
#' of the functions implemented in the present package.
#'
#' Each function is tested during the building of the package,
#' which means that results are guaranteed to match those of
#' the equivalent Matlab functions to at least 8 digits.
#'
#' A significant difference from the Matlab case is in the inspection
#' of the dimensions of arguments. The Matlab library has rules 
#' for expanding some arguments to match others. For example,
#' if Practical Salinity is a matrix and pressure is a single value,
#' then that single pressure is used throughout a calculation of 
#' Absolute Salinity. This convenience is only partly mimicked in the
#' present package.  Since the underlying C code works on vectors, 
#' the R functions in \code{gsw} start by transforming the arguments accordingly.
#' This involves using \code{\link{rep}} on each argument to get something
#' with length matching the first argument, and, after the computation
#' is complete, converting the return value into a matrix, if the first
#' argument was a matrix. There are some exceptions to this, however.
#' For example, \code{\link{gsw_SA_from_SP}} and similar functions
#' can handle the case in which the \code{SA} argument is a matrix and
#' \code{longitude} and \code{latitude} are vectors sized to match. 
#' This can be handy with gridded datasets. However, the careful
#' analyst will probably prefer to avoid this and other conveniences,
#' supplying properly-matched arguments from the outset.
#'
#' @docType package
#' @name gsw
NULL

#' Global SA lookup file
#'
#' @description
#' This dataset is not intended for users, but rather for internal use
#' within the \code{gsw} package. The dataset stores the 1.4M lookup
#' table defined in the 8.3M file \code{src/gsw_saar_data.c} in the C
#' library. (The .c file exceeds CRAN limitations on size.)
#'
#' @details
#' The data are designed to replace C elements defined as below
#' in \code{src/gsw_saar_data.c}:
#' \preformatted{
#'     static int	gsw_nx=91, gsw_ny=45, gsw_nz=45;
#'     static double	longs_ref[91];
#'     static double	lats_ref[45];
#'     static double	p_ref[45];
#'     static double	ndepth_ref[4095];
#'     static double	saar_ref[184275];
#'     static double	delta_sa_ref[184275];
#' }
#' 
#' R storage is in a list named \code{saar}, with elements named
#' as in the C code, i.e. \code{gsw_nx} etc.
#'
#' C storage for these variables is allocated as needed,
#' and the data are inserted, when \code{gsw} is launched.
#' Thus, the existing C library code "knows" about the data
#' as local storage, which keeps alterations to the C library to 
#' a minimum.
#'
#' The code used to create the RDA file (using the Fortran data
#' file, version 3.0.3) is given below.
#' \preformatted{
#'     gsw_nx <- 91
#'     gsw_ny <- 45
#'     gsw_nz <- 45
#'     f <- file("~/src/gsw_fortran_v3_03/gsw_data_v3_0.dat", "r")
#'     longs_ref <- scan(f, double(), n=gsw_nx)
#'     lats_ref <- scan(f, double(), n=gsw_ny)
#'     p_ref <- scan(f, double(), n=gsw_nz)
#'     ndepth_ref <- scan(f, double(), n=gsw_nx*gsw_ny)
#'     saar_ref <- scan(f, double(), n=gsw_nx*gsw_ny*gsw_nz)
#'     delta_sa_ref <- scan(f, double(), n=gsw_nx*gsw_ny*gsw_nz)
#'     saar <- list(gsw_nx=gsw_nx, gsw_ny=gsw_ny, gsw_nz=gsw_nz,
#'                  longs_ref=longs_ref, lats_ref=lats_ref, p_ref=p_ref, ndepth_ref=ndepth_ref,
#'                  saar_ref=saar_ref, delta_sa_ref=delta_sa_ref)
#'     save(saar, file="saar.rda")
#'     tools::resaveRdaFiles("saar.rda")
#'     close(f)
#'}
#'
#' @docType data
#' @name saar
NULL


## PART 2: utility functions

#' Reshape list elements to match the shape of the first element.
#'
#' This is mainly used within gsw, to ensure that arguments sent
#' to the C functions are of equal length.  This is a convenience, 
#' for processing data that often have this condition. For example, a
#' CTD profile is likely to have many values for SP, t, and p,
#' but just a single value for each of longitude and latitude.
#' It is important to call argfix() to handle such cases, because
#' otherwise the underlying C code will be looking past the end of
#' the vectors storing longitude and latitude, which can yield odd
#' results or even segmentation faults.
#'
#' @param list A list of elements, typically arguments that will be used in GSW functions.
#' @return A list with all elements of same shape (length or dimension).
argfix <- function(list)
{
    n <- length(list)
    if (n > 0) {
        length1 <- length(list[[1]])
        for (i in 2:n) {
            if (length(list[[i]]) != length1) {
                list[[i]] <- rep(list[[i]], length.out=length1)
            }
        }
        if (is.matrix(list[[1]])) {
            for (i in 2:n) {
                dim(list[[i]]) <- dim(list[[1]])
            }
        }
    }
    list
}



## PART 3: gsw (Gibbs SeaWater) functions, in alphabetical order (ignoring case)

#' Adiabatic lapse rate from Conservative Temperature
#'
#' Note that the unit is K/Pa, i.e. 1e-4 times K/dbar.
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return adiabatic lapse rate (note unconventional unit) [ K/Pa ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' lr <- gsw_adiabatic_lapse_rate_from_CT(SA, CT, p)
#' expect_equal(lr, 1e-7*c(0.240199646230069, 0.238457486976761, 0.203635157319712,
#'                         0.119829566859790, 0.100052760967308, 0.087773070307283))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_adiabatic_lapse_rate_from_CT.html}
gsw_adiabatic_lapse_rate_from_CT <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_adiabatic_lapse_rate_from_CT",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}
                                        
#' Thermal expansion coefficient with respect to Conservative Temperature. (75-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return thermal expansion coefficient with respect to Conservative Temperature [ 1/K ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' alpha <- gsw_alpha(SA,CT,p)
#' expect_equal(alpha, 1e-3 * c( 0.324464211877393, 0.322610094680523, 0.281335030247435,
#'                              0.173529986885424, 0.146898108553385, 0.130265123640082))
#' @seealso The salinity analogue to this is \code{\link{gsw_beta}}; other related functions include \code{\link{gsw_beta_const_t_exact}}, \code{\link{gsw_alpha_wrt_t_exact}} and \code{\link{gsw_alpha_on_beta}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_alpha.html}
gsw_alpha <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_alpha",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Thermal expansion coefficient over haline contraction coefficient. (75-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return ratio of thermal expansion coefficient to haline contraction coefficient [ (g/kg)/K ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' alpha_on_beta <- gsw_alpha_on_beta(SA,CT,p)
#' expect_equal(alpha_on_beta, c(0.452468543022009, 0.449601695030057, 0.387140203094424,
#'                               0.230778871228268, 0.193747796234162, 0.170946048860385))
#' @seealso This yields the ratio of the return values from \code{\link{gsw_alpha}} and \code{\link{gsw_beta}}, to within computational precision.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_alpha_on_beta.html}
gsw_alpha_on_beta <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_alpha_on_beta",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Thermal expansion coefficient with respect to in-situ temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90)  [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return thermal expansion coefficient with respect to in-situ temperature [ 1/K ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <-  c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' alpha_wrt_t_exact <- gsw_alpha_wrt_t_exact(SA,t,p)
#' expect_equal(alpha_wrt_t_exact, 1e-3*c(0.325601747227247, 0.323448083851267, 0.281413883319329,
#'                                        0.172825692975230, 0.145569941503599, 0.128362986933288))
#' 
#' @seealso \code{\link{gsw_alpha}}, \code{\link{gsw_beta}} and \code{\link{gsw_alpha_on_beta}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_alpha_wrt_t_exact.html}
gsw_alpha_wrt_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_alpha_wrt_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Saline contraction coefficient at constant Conservative Temperature. (75-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return saline contraction coefficient at constant Conservative Temperature [ kg/g ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' beta <- gsw_beta(SA,CT,p)
#' expect_equal(beta, 1e-3*c(0.717521909550091, 0.717657376442386, 0.726169785748549,
#'                           0.750420924314564, 0.754903052075032, 0.756841573481865))
#' 
#' @seealso
#' The temperature analogue to this is \code{\link{gsw_alpha}}; other related functions
#' include \code{\link{gsw_alpha_wrt_t_exact}} and \code{\link{gsw_alpha_on_beta}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_beta.html}
gsw_beta <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_beta",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Saline contraction coefficient at constant in-situ temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return saline contraction coefficient at constant in-situ temperature [ kg/g ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <- c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' b <- gsw_beta_const_t_exact(SA, t, p)
#' expect_equal(b, 1e-3*c(0.731120837010429, 0.731071779078011, 0.736019128913071,
#'                        0.753810501711847, 0.757259405338257, 0.758649268096996))
#' @seealso
#' A related function is \code{\link{gsw_beta}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_beta_const_t_exact.html}
gsw_beta_const_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_beta_const_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Cabbeling coefficient (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return cabbeling coefficient with respect to Conservative Temperature [ 1/(K^2) ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' cabbeling <- gsw_cabbeling(SA,CT,p)
#' expect_equal(cabbeling, 1e-4*c(0.086645721047423, 0.086837829466794, 0.092525582052438,
#'                                0.108884336975401, 0.112971197222338, 0.115483896148927))
#' 
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_cabbeling.html}
gsw_cabbeling <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_cabbeling",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}
#' Electrical conductivity from Practical Salinity
#'
#' Note: the return value is not conductivity ratio, but rather
#' conductivity itself, in mS/cm.  To convert to conductivity ratio,
#' divide by 42.9140 (the value of conductivity at S=35, T68=15, and
#' p=0).
#' 
#' @param SP Practical Salinity (PSS-78) [ unitless ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return electrical conductivity [ mS/cm ]
#' @examples 
#' library(testthat)
#' SP <- c(34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600)
#' t <- c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' C <- gsw_C_from_SP(SP, t, p)
#' expect_equal(C, c(56.412599581571186, 56.316185602699953, 50.670369333973944,
#'                   38.134518936104350, 35.056577637635257, 32.986550607990118))
#' @seealso \code{\link{gsw_SP_from_C}} does the reverse.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_C_from_SP.html}
gsw_C_from_SP <- function(SP, t, p)
{
    l <- argfix(list(SP=SP, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_C_from_SP",
               SP=as.double(l$SP), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Isobaric heat capacity
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return heat capacity [ J/(kg*K) ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <-  c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' cp_t_exact <- gsw_cp_t_exact(SA, t, p)
#' expect_equal(cp_t_exact, 1e3*c(4.002888003958537, 4.000980283927373, 3.995546468894633,
#'                                3.985076769021370, 3.973593843482723, 3.960184084786622))
#' 
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_cp_t_exact.html}
gsw_cp_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_cp_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}


#' Conservative temperature freezing point
#'
#' Note: as of 2014-12-23, this corresponds to the Matlab function
#' called \code{gsw_t_freezing_poly}. (The confusion arises from a
#' mismatch in release version between the Matlab and C libraries.)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param saturation_fraction saturation fraction of dissolved air in seawater
#' @return Conservative Temperature at freezing of seawater [ deg C ]. That is, the freezing temperature expressed in terms of Conservative Temperature (ITS-90). 
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' saturation_fraction <- 1
#' CT_freezing <- gsw_CT_freezing(SA, p, saturation_fraction)
#' expect_equal(CT_freezing, c(-1.899683776424096, -1.940791867869104, -2.006240664432488,
#'                             -2.092357761318778, -2.359300831770506, -2.677162675412748))
#' @seealso \code{\link{gsw_t_freezing}} is the analogue for in-situ temperature.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_CT_freezing.html}
gsw_CT_freezing <- function(SA, p, saturation_fraction=1)
{
    l <- argfix(list(SA=SA, p=p, saturation_fraction=saturation_fraction))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_CT_freezing",
               SA=as.double(l$SA), p=as.double(l$p), saturation_fraction=as.double(l$saturation_fraction),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Conservative Temperature from potential temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param pt potential temperature (ITS-90) [ deg C ]
#' @return Conservative Temperature [ deg C ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' pt <- c(28.7832, 28.4209, 22.7850, 10.2305,  6.8292,  4.3245)
#' CT <- gsw_CT_from_pt(SA, pt)
#' expect_equal(CT, c(28.809923015982083, 28.439144260767169, 22.786246608464264,
#'                    10.226165605435785, 6.827183417643142,  4.323565182322069))
#' @seealso
#' \code{\link{gsw_CT_from_t}} calculates Conservative Temperature from in-situ temperature
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_CT_from_pt.html}
gsw_CT_from_pt <- function(SA, pt)
{
    l <- argfix(list(SA=SA, pt=pt))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_CT_from_pt",
               SA=as.double(l$SA), pt=as.double(l$pt),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Convert from temperature to conservative temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Conservative Temperature [ deg C ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <-  c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' CT <- gsw_CT_from_t(SA, t, p)
#' expect_equal(CT, c(28.809919826700281, 28.439227816091140, 22.786176893078498,
#'                    10.226189266620782, 6.827213633479988, 4.323575748610455))
#' 
#' @seealso \code{\link{gsw_t_from_CT}} does the reverse
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_CT_from_t.html}
gsw_CT_from_t <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_CT_from_t",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Absolute Salinity Anomaly from Practical Salinity
#' 
#' @param SP Practical Salinity  (PSS-78) [ unitless ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
#' @param latitude latitude in decimal degrees north [ -90 ... +90 ]
#' @return deltaSA Absolute Salinity Anomaly  [ g/kg ]
#' @examples 
#' library(testthat)
#' SP =   c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' p =    c(     10,      50,     125,     250,     600,    1000)
#' lat =  c(      4,       4,       4,       4,       4,       4)
#' long = c(    188,     188,     188,     188,     188,     188)
#' deltaSA = gsw_deltaSA_from_SP(SP,p,long,lat)
#' expect_equal(deltaSA, c(0.000167203365230, 0.000268836122231, 0.000665803155705,
#'                         0.002706154619403, 0.005652977406832,  0.009444734661606))
#' @seealso
#' \code{\link{gsw_SA_from_SP}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_deltaSA_from_SP.html}
gsw_deltaSA_from_SP <- function(SP, p, longitude, latitude)
{
    if (missing(SP)) stop("must supply SP")
    if (missing(p)) stop("must supply p")
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    l <- argfix(list(SP=SP, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_deltaSA_from_SP",
               SP=as.double(l$SP), p=as.double(l$p),
               longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Dynamic enthalpy of seawater (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return dynamic enthalpy [ J/kg ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <-c( 28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' de <- gsw_dynamic_enthalpy(SA, CT, p)
#' expect_equal(de, 1e3*c(0.097864698087770, 0.489161476686235, 1.220512192086506,
#'                        2.433731199531144, 5.833880057399701, 9.711443860944032))
#' @seealso
#' \code{\link{gsw_enthalpy}} and \code{\link{gsw_enthalpy_t_exact}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_enthalpy.html}
gsw_dynamic_enthalpy <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_dynamic_enthalpy",
               SA=as.double(l$SA), t=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Specific enthalpy of seawater (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return specific enthalpy [ J/kg ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <-c( 28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' e <- gsw_enthalpy(SA, CT, p)
#' expect_equal(e, 1e5*c(1.151031813559086, 1.140146926828028, 0.921800138366058,
#'                       0.432553713026279, 0.330871609742468, 0.269706841603465))                  
#' @seealso
#' \code{\link{gsw_dynamic_enthalpy}} and \code{\link{gsw_enthalpy_t_exact}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_enthalpy.html}
gsw_enthalpy <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_enthalpy",
               SA=as.double(l$SA), t=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Specific enthalpy of seawater
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90)  [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return specific enthalpy [ J/kg ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <- c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' e <- gsw_enthalpy_t_exact(SA, t, p)
#' expect_equal(e, 1e5*c(1.151032604783763, 1.140148036012021, 0.921799209310966,
#'                       0.432553283808897, 0.330872159700175, 0.269705880448018))
#' @seealso
#' \code{\link{gsw_enthalpy}} and \code{\link{gsw_dynamic_enthalpy}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_enthalpy_t_exact.html}
gsw_enthalpy_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_enthalpy_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Specific entropy as a function of in-situ temperature and pressure
#'
#' Calculates specific entropy given Absolute Salinity, in-situ
#' temperature and pressure.
#'
#' The related function gsw_entropy_from_CT() is not provided
#' in the C library, although it is available in the (later-
#' versioned) Matlab library.
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return specific entropy [ J/(kg*K) ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <- c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' e <- gsw_entropy_from_t(SA, t, p)
#' expect_equal(e, 100*c(4.003894252787245, 3.954381784340642, 3.198664981986740,
#'                       1.467908815899072, 0.986473408657975, 0.627915087346090))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_entropy_from_t.html}
gsw_entropy_from_t <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_entropy_from_t",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Gravitational acceleration
#' 
#' @param latitude latitude in decimal degress north [ -90 ... +90 ]
#' @param p sea pressure [ dbar ]
#' @return gravitational acceleration [ m/s^2 ]
#' @examples
#' library(testthat)
#' lat <- c(-90, -60, -30, 0)
#' grav <- gsw_grav(lat)
#' expect_equal(grav, c(9.832186205884799, 9.819178859991149,
#'                      9.793249257048750, 9.780327000000000))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_grav.html}
gsw_grav <- function(latitude, p=0)
{
    l <- argfix(list(latitude=latitude, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_grav",
               latitude=as.double(l$latitude), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(latitude))
        dim(rval) <- dim(latitude)
    rval
}

#' Specific internal energy of seawater (75-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return specific internal energy [ J/kg ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' e <- gsw_internal_energy(SA, CT, p)
#' expect_equal(e, 1e5*c(1.148091576956162, 1.134013145527675, 0.909571141498779,
#'                       0.408593072177020, 0.273985276460357, 0.175019409258405))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_internal_energy.html}
gsw_internal_energy <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_internal_energy",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Ratio of vert. gradient of pot. density to vert grad of locally-referenced pot density
#'
#' Note that the C library had to be patched to get this working; a new
#' version of the library will address the bug directly.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @param p_ref reference pressure [ dbar ]
#' @return list containing IPV_vs_fNsquared_ratio [ unitless ] and mid-point pressure p_mid [ dbar ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' p_ref <- 0
#' r <- gsw_IPV_vs_fNsquared_ratio(SA, CT, p, p_ref)
#' expect_equal(r$IPV_vs_fNsquared_ratio, c(0.999742244888022, 0.996939883468178, 0.986141997098021,
#'                                          0.931595598713477, 0.861224354872028))
#' expect_equal(r$p_mid, c(30, 87.5, 187.5, 425, 800))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_IPV_vs_fNsquared_ratio.html}
gsw_IPV_vs_fNsquared_ratio <- function(SA, CT, p, p_ref=0)
{
    l <- argfix(list(SA=SA, CT=CT, p=p, p_ref=p_ref))
    n <- length(l[[1]])
    ## note: only use p_ref[1] since the C-library code says it must be constant
    r <- .C("wrap_gsw_IPV_vs_fNsquared_ratio",
            SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p), p_ref=as.double(l$p_ref[1]),
            n=as.integer(n),
            ratio=double(n-1), p_mid=double(n-1), NAOK=TRUE, package="gsw")
    if (is.matrix(SA))
        stop("gsw_IPV_vs_fNsquared_ratio() cannot handle matrix SA")
    list(IPV_vs_fNsquared_ratio=r$ratio, p_mid=r$p_mid)
}

#' Isentropic compressibility of seawater (75-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return isentropic compressibility [ 1/Pa ] (not 1/dbar)
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <-c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' kappa <- gsw_kappa(SA, CT, p)
#' expect_equal(kappa, 1e-9*c(0.411343648791300, 0.411105416128094, 0.416566236026610,
#'                            0.435588650838751, 0.438782500588955, 0.439842289994702))
#' @seealso \code{\link{gsw_kappa_t_exact}} is an analogue in terms of in-situ temperature
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_kappa.html}
gsw_kappa <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_kappa",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Isentropic compressibility of seawater
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return isentropic compressibility [ 1/Pa ] (not 1/dbar)
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <-c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' kappa <- gsw_kappa(SA, CT, p)
#' expect_equal(kappa, 1e-9*c(0.411343648791300, 0.411105416128094, 0.416566236026610,
#'                            0.435588650838751, 0.438782500588955, 0.439842289994702))
#' @seealso \code{\link{gsw_kappa}} is an analogue in terms of Conservative Temperature
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_kappa_t_exact.html}
gsw_kappa_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_kappa_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Latent heat of evaporation
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return latent heat of evaporation [ J/kg ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' lh  <- gsw_latentheat_evap_CT(SA, CT)
#' expect_equal(lh, 1e6*c(2.429947107462561, 2.430774073049213, 2.444220372158452,
#'                        2.474127109232524, 2.482151446148560, 2.488052297193594))
#' @seealso \code{\link{gsw_latentheat_evap_t}} is an analogue in terms of in-situ temperature. For melting, see \code{\link{gsw_latentheat_melting}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_latentheat_evap_CT.html}
gsw_latentheat_evap_CT <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_latentheat_evap_CT",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Latent heat of evaporation
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @return latent heat of evaporation [ J/kg ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <- c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' lh = gsw_latentheat_evap_t(SA, t)
#' expect_equal(lh, 1e6*c(2.429882982734836, 2.430730236218543, 2.444217294049004,
#'                        2.474137411322517, 2.482156276375029, 2.488054617630297))
#' @seealso \code{\link{gsw_latentheat_evap_CT}} is an analogue in terms of Conservative Temperature. For melting, see \code{\link{gsw_latentheat_melting}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_latentheat_evap_t.html}
gsw_latentheat_evap_t <- function(SA, t)
{
    l <- argfix(list(SA=SA, t=t))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_latentheat_evap_t",
               SA=as.double(l$SA), t=as.double(l$t),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Latent heat of melting
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @return latent heat of freezing [ J/kg ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' lh <- gsw_latentheat_melting(SA, p)
#' expect_equal(lh, 1e5*c(3.299495187300804, 3.298611954422526, 3.297124383647952,
#'                        3.294972884747496, 3.288480015369891, 3.280715953443947),
#'              scale=1e5, tolerance=1e-6)
#' # FIXME: why don't scale=NULL and tolerance=.Machine$double.eps^0.5 work?
#' @seealso \code{\link{gsw_latentheat_evap_CT}} and \code{\link{gsw_latentheat_evap_t}} are analogues for evaporation.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_latentheat_melting.html}
gsw_latentheat_melting <- function(SA, p)
{
    l <- argfix(list(SA=SA, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_latentheat_melting",
               SA=as.double(l$SA), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Calculate Brunt Vaisala Frequency squared
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return list containing N2 [ 1/s^ ] and mid-point pressure p_mid [ dbar ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' latitude <- 4
#' r <- gsw_Nsquared(SA, CT, p, latitude)
#' N2 <- r$N2
#' p_mid <- r$p_mid
#' expect_equal(N2, 1e-3*c(0.060843209693499, 0.235723066151305, 0.216599928330380,
#'                         0.012941204313372, 0.008434782795209))
#' expect_equal(p_mid, c(30, 87.5, 187.5, 425, 800))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html}
gsw_Nsquared <- function(SA, CT, p, latitude=0)
{
    l <- argfix(list(SA=SA, CT=CT, p=p, latitude=latitude))
    n <- length(l[[1]])
    r <- .C("wrap_gsw_Nsquared",
            SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p), latitude=as.double(l$latitude),
            n=as.integer(n), n2=double(n-1), p_mid=double(n-1), NAOK=TRUE, package="gsw")
    if (is.matrix(SA))
        stop("gsw_Nsquared() cannot handle matrix SA")
    list(N2=r$n2, p_mid=r$p_mid)
}

#' Pressure from height (75-term equation)
#' 
#' @param z height, zero at surface (but note last 2 args) and positive upwards [ m ]
#' @param latitude latitude in decimal degrees north [ -90 ... +90 ]
#' @param geo_strf_dyn_height dynamic height anomaly [ m^2/s^2 ]
#' @param sea_surface_geopotential geopotential at zero sea pressure [ m^2/s^2 ]
#' @return sea pressure [ dbar ]
#' @examples
#' library(testthat)
#' z <- -c(10, 50, 125, 250, 600, 1000)
#' latitude <- 4
#' p <- gsw_p_from_z(z, latitude)
#' expect_equal(p, 1e3*c(0.010055726724518, 0.050283543374874, 0.125731858435610,
#'                       0.251540299593468, 0.604210012340727, 1.007990337692001))
#' @seealso
#' This is (almost) the reverse of \code{\link{gsw_z_from_p}}, apart from the last two arguments.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_p_from_z.html}
gsw_p_from_z <- function(z, latitude, geo_strf_dyn_height=0, sea_surface_geopotential=0)
{
    l <- argfix(list(z=z, latitude=latitude,
                     geo_strf_dyn_height=geo_strf_dyn_height,
                     sea_surface_geopotential=sea_surface_geopotential))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_p_from_z",
               z=as.double(l$z), latitude=as.double(l$latitude),
               geo_strf_dyn_height=as.double(l$geo_strf_dyn_height),
               sea_surface_geopotential=as.double(l$sea_surface_geopotential),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(z))
        dim(rval) <- dim(z)
    rval
}

#' Potential density
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @param p_ref reference pressure [ dbar ]
#' @return potential density [ kg/m^3 ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <- c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' p_ref  <- 0
#' prho <- gsw_pot_rho_t_exact(SA,t,p,p_ref)
#' expect_equal(prho, 1e3*c(1.021798145811089, 1.022052484416980, 1.023893583651958,
#'                          1.026667621124443, 1.027107230868492, 1.027409631264134)) 
#' @seealso
#' \code{\link{gsw_rho}} and \code{\link{gsw_rho_t_exact}} compute density; \code{\link{gsw_sigma0}} and related functions compute potential density at particular pressures.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_pot_rho_t_exact.html}
gsw_pot_rho_t_exact <- function(SA, t, p, p_ref)
{
    l <- argfix(list(SA=SA, t=t, p=p, p_ref=p_ref))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_pot_rho_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p), pref=as.double(l$p_ref),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential temperature referenced to the surface
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return potential temperature [ deg C ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <- c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' pt0 <- gsw_pt0_from_t(SA, t, p)
#' expect_equal(pt0, c(28.783196819670632, 28.420983342398962, 22.784930399117108,
#'                     10.230523661095731, 6.829230224409661, 4.324510571845719))
#' @seealso \code{\link{gsw_pt_from_CT}} and \code{\link{gsw_pt_from_t}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_pt0_from_t.html}
gsw_pt0_from_t <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_pt0_from_t",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential temperature from Conservative Temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential temperature [ deg C ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' pt <- gsw_pt_from_CT(SA, CT)
#' expect_equal(pt, c(28.783177048624573, 28.420955597191984, 22.784953468087107,
#'                    10.230534394434429, 6.829216587061605, 4.324534835990236))
#' @seealso \code{\link{gsw_pt0_from_t}} for the surface case and and \code{\link{gsw_pt_from_t}} for the general case.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_pt_from_CT.html}
gsw_pt_from_CT <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_pt_from_CT",
               SA=as.double(l$SA), t=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential temperature from in-situ temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90)  [ deg C ]
#' @param p sea pressure [ dbar ]
#' @param p_ref reference pressure [ dbar ]
#' @return potential temperature [ deg C ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <- c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' p_ref <- 0
#' pt <- gsw_pt_from_t(SA, t, p, p_ref)
#' expect_equal(pt, c(28.783196819670632, 28.420983342398962, 22.784930399117108,
#'                    10.230523661095731, 6.829230224409661, 4.324510571845719)) 
#' @seealso \code{\link{gsw_pt_from_CT}} is the analogue for Conservative Temperature
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_pt_from_t.html}
gsw_pt_from_t <- function(SA, t, p, p_ref=0)
{
    l <- argfix(list(SA=SA, t=t, p=p, p_ref=p_ref))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_pt_from_t",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p), p_ref=as.double(l$p_ref),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' In-situ density (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return in-situ density [ kg/m^3 ]
#' @seealso \code{\link{gsw_rho_t_exact}} is similar to this, but using in-situ temperature. SA and CT may be computed from UNESCO quantities using \code{\link{gsw_SA_from_SP}} and \code{\link{gsw_CT_from_t}}. For potential density anomalies, use \code{\link{gsw_sigma0}} and related functions.
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' rho <- gsw_rho(SA,CT,p)
#' expect_equal(rho, 1e3*c(1.021839935738108, 1.022262457966867, 1.024427195413316,
#'                         1.027790152759127, 1.029837779000189, 1.032002453224572))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_rho.html}
gsw_rho <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_rho",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' SA, CT and p partial derivatives of density (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return list containing drho_dSA [ kg^2/(g m^3) ], drho_dCT [ kg/(K m^3) ] and drho_dp [ kg/(Pa m^3) ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' r <- gsw_rho_first_derivatives(SA, CT, p)
#' expect_equal(r$drho_dSA, c(0.733153791778356, 0.733624109867480, 0.743950957375504,
#'                            0.771357282286743, 0.777581141431288, 0.781278296628328))
#' expect_equal(r$drho_dCT, c(-0.331729027977015, -0.329838643311336, -0.288013324730644,
#'                            -0.178012962919839, -0.150654632545556, -0.133556437868984))
#' expect_equal(r$drho_dp, 1e-6*c(0.420302360738476, 0.420251070273888, 0.426773054953941,
#'                                0.447763615252861, 0.452011501791479, 0.454118117103094))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_rho_first_derivatives.html}
gsw_rho_first_derivatives <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_rho_first_derivatives",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, drho_dSA=double(n), drho_dCT=double(n), drho_dp=double(n),
               NAOK=TRUE, package="gsw")
    if (is.matrix(SA))
        stop("gsw_rho_first_derivatives() cannot handle matrix SA")
    list(drho_dSA=rval$drho_dSA, drho_dCT=rval$drho_dCT, drho_dp=rval$drho_dp)
}

#' In-situ density
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return in-situ density [ kg/m^3 ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <-  c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' rho <- gsw_rho_t_exact(SA, t, p)
#' expect_equal(rho, 1e3*c(1.021840173185531, 1.022262689926782, 1.024427715941676,
#'                         1.027790201811623, 1.029837714725961, 1.032002404116447))
#' @seealso \code{\link{gsw_rho}} is similar but uses SA and CT; SA may be computed from UNESCO quantities using \code{\link{gsw_SA_from_SP}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_rho_t_exact.html}
gsw_rho_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_rho_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Convert from density to absolute salinity
#'
#' @param rho seawater density [ kg/m^3 ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Absolute Salinity [ g/kg ]
#' @examples
#' library(testthat)
#' rho <- c(1021.8482, 1022.2647, 1024.4207, 1027.7841, 1029.8287, 1031.9916)
#' CT <-c(    28.7856,   28.4329,   22.8103,   10.2600,    6.8863,    4.4036)
#' p <- c(         10,        50,       125,       250,       600,      1000)
#' SA <- gsw_SA_from_rho(rho, CT, p)
#' expect_equal(SA, c(34.712080120418108, 34.891723808488869, 35.026202257609505,
#'                    34.847160842234572, 34.736398269039945, 34.732228881079742))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_rho.html}
gsw_SA_from_rho <- function(rho, CT, p)
{
    l <- argfix(list(rho=rho, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SA_from_rho",
               SA=as.double(l$rho), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(rho))
        dim(rval) <- dim(rho)
    rval
}

#' Convert from practical salinity to absolute salinity
#'
#' Calculate Absolute Salinity from Practical Salinity, pressure,
#' longitude, and latitude.
#'
#' If SP is a matrix and if its dimensions correspond to the
#' lengths of longitude and latitude, then the latter are
#' converted to analogous matrices with \code{\link{expand.grid}}.
#' 
#' @param SP Practical Salinity (PSS-78) [ unitless ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return Absolute Salinity [ g/kg ]
#' @examples
#' library(testthat)
#' SP <- c(34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' lat <- c(     4,       4,       4,       4,       4,       4)
#' long <- c(  188,     188,     188,     188,     188,     188)
#' SA <- gsw_SA_from_SP(SP, p, long, lat)
#' expect_equal(SA, c(34.711778344814114, 34.891522618230098, 35.025544862476920,
#'                    34.847229026189588, 34.736628474576051, 34.732363065590846))
#' @seealso \code{\link{gsw_SP_from_SA}} does the reverse.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_SP.html}
gsw_SA_from_SP <- function(SP, p, longitude, latitude)
{
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    ## check for special case that SP is a matrix defined on lon and lat
    if (is.matrix(SP)) {
        dim <- dim(SP)
        if (length(longitude) == dim[1] && length(latitude) == dim[2]) {
            ll <- expand.grid(longitude=as.vector(longitude), latitude=as.vector(latitude))
            longitude <- ll$longitude
            latitude <- ll$latitude
        }
    }
    l <- argfix(list(SP=SP, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SA_from_SP",
               SP=as.double(l$SP), p=as.double(l$p),
               longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Absolute Salinity from Preformed Salinity
#'
#' Calculate Absolute Salinity from Preformed Salinity, pressure,
#' longitude, and latitude.
#'
#' If Sstar is a matrix and if its dimensions correspond to the
#' lengths of longitude and latitude, then the latter are
#' converted to analogous matrices with \code{\link{expand.grid}}.
#' 
#' @param Sstar Preformed Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return Absolute Salinity [ g/kg ]
#' @examples
#' library(testthat)
#' SP <- c(34.7115, 34.8912, 35.0247, 34.8436, 34.7291, 34.7197)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' lat <- c(     4,       4,       4,       4,       4,       4)
#' long <- c(  188,     188,     188,     188,     188,     188)
#' SA <- gsw_SA_from_Sstar(SP, p, long, lat)
#' expect_equal(SA, c(34.711724663585905, 34.891561223296009, 35.025594598699882,
#'                    34.847235885385913, 34.736694493054166, 34.732387111902753))
#' @seealso \code{\link{gsw_Sstar_from_SA}} does the reverse.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_Sstar.html}
gsw_SA_from_Sstar <- function(Sstar, p, longitude, latitude)
{
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    ## check for special case that Sstar is a matrix defined on lon and lat
    if (is.matrix(Sstar)) {
        dim <- dim(Sstar)
        if (length(longitude) == dim[1] && length(latitude) == dim[2]) {
            ll <- expand.grid(longitude=as.vector(longitude), latitude=as.vector(latitude))
            longitude <- ll$longitude
            latitude <- ll$latitude
        }
    }
    l <- argfix(list(Sstar=Sstar, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SA_from_Sstar",
               Sstar=as.double(l$Sstar), p=as.double(l$p),
               longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(Sstar))
        dim(rval) <- dim(Sstar)
    rval
}

#' Potential density anomaly referenced to 0 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 0 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly [ kg/m^3 ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' sigma0 <- gsw_sigma0(SA,CT)
#' expect_equal(sigma0, c(21.797900819337656, 22.052215404397316, 23.892985307893923,
#'                        26.667608665972011, 27.107380455119710, 27.409748977090885))
#' @seealso Use \code{\link{gsw_sigma1}} for 1000 dbar pressure, \code{\link{gsw_sigma2}} for 2000 dbar, \code{\link{gsw_sigma3}} for 3000 dbar, or \code{\link{gsw_sigma4}} for 4000 dbar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma0.html}
gsw_sigma0 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma0",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential density anomaly referenced to 1000 dbar
#'
#' This uses the 75-term density equation, and returns
#' potential density referenced to a pressure of 1000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly [ kg/m^3 ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' sigma1 <- gsw_sigma1(SA,CT)
#' expect_equal(sigma1, c(25.955618850310202, 26.213131422420247, 28.125423775188438,
#'                        31.120360038882382, 31.637724222733368, 32.002453224572037))
#' @seealso Use \code{\link{gsw_sigma0}} for 0 dbar pressure, \code{\link{gsw_sigma2}} for 2000 dbar, \code{\link{gsw_sigma3}} for 3000 dbar, or \code{\link{gsw_sigma4}} for 4000 dbar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma1.html}
gsw_sigma1 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma1",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential density anomaly referenced to 2000 dbar
#'
#' This uses the 75-term density equation, and returns
#' potential density referenced to a pressure of 2000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly [ kg/m^3 ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' sigma2 <- gsw_sigma2(SA,CT)
#' expect_equal(sigma2, c(30.023152223799116, 30.283783336283477, 32.265556840289719,
#'                        35.474550881051073, 36.067289438047737, 36.492606494879510)) 
#' @seealso Use \code{\link{gsw_sigma0}} for 0 dbar pressure, \code{\link{gsw_sigma1}} for 1000 dbar, \code{\link{gsw_sigma3}} for 3000 dbar, or \code{\link{gsw_sigma4}} for 4000 dbar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma2.html}
gsw_sigma2 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma2",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential density anomaly referenced to 3000 dbar
#'
#' This uses the 75-term density equation, and returns
#' potential density referenced to a pressure of 3000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly with reference pressure 3000 dbar [ kg/m^3 ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' sigma3 <- gsw_sigma3(SA,CT)
#' expect_equal(sigma3, c(34.003747849903675, 34.267409891564057, 36.316415829697917,
#'                        39.732367693977039, 40.397934186745033, 40.881795690566832))
#' @seealso Use \code{\link{gsw_sigma0}} for 0 dbar pressure, \code{\link{gsw_sigma1}} for 1000 dbar, \code{\link{gsw_sigma2}} for 2000 dbar, or \code{\link{gsw_sigma4}} for 4000 dbar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma3.html}
gsw_sigma3 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma3",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Potential density anomaly referenced to 4000 dbar
#'
#' This uses the 75-term density equation, and returns
#' potential density referenced to a pressure of 4000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly with reference pressure 4000 dbar [ kg/m^3 ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' sigma4 <- gsw_sigma4(SA,CT)
#' expect_equal(sigma4, c(37.900374609834898, 38.166979617032439, 40.280876075282549,
#'                        43.896091033421953, 44.631677245327637, 45.171817312020039))
#' @seealso Use \code{\link{gsw_sigma0}} for 0 dbar pressure, \code{\link{gsw_sigma1}} for 1000 dbar, \code{\link{gsw_sigma2}} for 2000 dbar, or \code{\link{gsw_sigma3}} for 3000 dbar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma4.html}
gsw_sigma4 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma4",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Sound speed (75-term equation)
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return sound speed [ m/s ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' sound_speed <- gsw_sound_speed(SA,CT,p)
#' expect_equal(sound_speed, 1e3*c(1.542426412426373, 1.542558891663385, 1.530801535436184,
#'                                 1.494551099295314, 1.487622786765276, 1.484271672296205))
#' @seealso \code{\link{gsw_sound_speed_t_exact}} for a precise formula using in-situ temperature
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sound_speed.html}
gsw_sound_speed<- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sound_speed",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Sound speed
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return sound speed [ m/s ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' sound_speed <- gsw_sound_speed_t_exact(SA,CT,p)
#' expect_equal(sound_speed, 1e3*c(1.542615803587414, 1.542703534065789, 1.530844979136360,
#'                                 1.494409996920661, 1.487377102518027, 1.483934609078705))
#' @seealso \code{\link{gsw_sound_speed}} for an approximate formula using CT
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sound_speed_t_exact.html}
gsw_sound_speed_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sound_speed_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Specific volume
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Specific volume (1/density)
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' specvol <- gsw_specvol(SA, CT, p)
#' expect_equal(specvol, 1e-3*c(0.978626852431313, 0.978222365701325, 0.976155264597929,
#'                              0.972961258011157, 0.971026719344908, 0.968989944622149))
#' @seealso With in-situ temperature, use \code{\link{gsw_specvol_t_exact}}; \code{\link{gsw_specvol_anom}} gives specific volume anomaly.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_specvol.html}
gsw_specvol  <- function(SA, CT, p)
{
    1 / gsw_rho(SA, CT, p)
}

#' Specific volume anomaly [standard] (75-term equation)
#'
#' Note that the TEOS function named \code{specific_volume_anomaly} is not
#' provided in the C library, so it is not provided in R, either.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' a <- gsw_specvol_anom_standard(SA, CT, p)
#' expect_equal(a, 1e-5*c(0.601051894897400, 0.578609769250563, 0.405600538950092,
#'                        0.142190453761838, 0.104335535578967, 0.076383389577725))
#' @return Specific volume anomaly [ m^3/kg ]
#' @seealso Specific volume itself is given by \code{\link{gsw_specvol}} and \code{\link{gsw_specvol_t_exact}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_specvol_anom_standard.html}
gsw_specvol_anom_standard <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_specvol_anom_standard", # FIXME: why the "standard" in name?
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Specific volume
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90)  [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Specific volume [ m^3/kg ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' t <- c( 28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' v <- gsw_specvol_t_exact(SA, t, p)
#' expect_equal(v, 1e-3 * c(0.978626625025472, 0.978222143734527, 0.976154768597586,
#'                          0.972961211575438, 0.971026779948624, 0.968989990731808))
#' @seealso With Conservative Temperature, use \code{\link{gsw_specvol}}; \code{\link{gsw_specvol_anom}} gives specific volume anomaly.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_specvol_t_exact.html}
gsw_specvol_t_exact  <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_specvol_t_exact",
               SA=as.double(l$SA), CT=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Convert from electrical conductivity to practical salinity
#' 
#' @param C conductivity [ mS/cm ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Practical Salinity (PSS-78) [ unitless ]
#' @examples 
#' library(testthat)
#' C <- c(34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600)
#' t <- c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p <- c(     10,      50,     125,     250,     600,    1000)
#' SP <- gsw_SP_from_C(C,t,p)
#' expect_equal(SP,  c(20.009869599086951, 20.265511864874270, 22.981513062527689,
#'                     31.204503263727982, 34.032315787432829, 36.400308494388170))
#' @seealso \code{\link{gsw_C_from_SP}} does the reverse.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_C.html}
gsw_SP_from_C <- function(C, t, p)
{
    l <- argfix(list(C=C, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SP_from_C",
               C=as.double(l$C), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(C))
        dim(rval) <- dim(C)
    rval
}

#' Convert from Absolute Salinity to Practical Salinity
#'
#' Calculate Practical Salinity from Absolute Salinity, pressure,
#' longitude, and latitude.
#'
#' If SP is a matrix and if its dimensions correspond to the
#' lengths of longitude and latitude, then the latter are
#' converted to analogous matrices with \code{\link{expand.grid}}.
#'
#' Note: unlike the corresponding Matlab function, this does not
#' return a flag indicating whether the location is in the ocean.
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return Practical Salinity (PSS-78) [ unitless ]
#' @examples 
#' library(testthat)
#' SA <-   c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' p <-    c(     10,      50,     125,     250,     600,    1000)
#' lat <-  c(      4,       4,       4,       4,       4,       4)
#' long <- c(    188,     188,     188,     188,     188,     188)
#' SP <- gsw_SP_from_SA(SA,p,long,lat)
#' expect_equal(SP, c(34.548721553448317, 34.727477488096639, 34.860554877708005,
#'                    34.680971112271791, 34.567971663653388, 34.560036751118204))
#' @seealso \code{\link{gsw_SA_from_SP}} does the reverse, while \code{\link{gsw_SP_from_SK}}, \code{\link{gsw_SP_from_SR}} and \code{\link{gsw_SP_from_Sstar}} are similar to this.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_SA.html}
gsw_SP_from_SA <- function(SA, p, longitude, latitude)
{
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    ## check for special case that SP is a matrix defined on lon and lat
    if (is.matrix(SA)) {
        dim <- dim(SA)
        if (length(longitude) == dim[1] && length(latitude) == dim[2]) {
            ll <- expand.grid(longitude=as.vector(longitude), latitude=as.vector(latitude))
            longitude <- ll$longitude
            latitude <- ll$latitude
        }
    }
    l <- argfix(list(SA=SA, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SP_from_SA",
               SA=as.double(l$SA), p=as.double(l$p), longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, SP=double(n), NAOK=TRUE, package="gsw")$SP
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Calculate Practical Salinity from Knudsen Salinity
#'
#' @param SK Knudsen Salinity [ parts per thousand, ppt ]
#' @return Practical Salinity (PSS-78) [ unitless ]
#' @examples 
#' library(testthat)
#' SK <- c(34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600)
#' SP <- gsw_SP_from_SK(SK)
#' expect_equal(SP, c(34.548342096952908, 34.727295637119113, 34.860409847645435,
#'                    34.680755706371187, 34.567658670360110, 34.559651800554022))
#' @seealso \code{\link{gsw_SP_from_SA}}, \code{\link{gsw_SP_from_SR}} and \code{\link{gsw_SP_from_Sstar}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_SK.html}
gsw_SP_from_SK <- function(SK)
{
    if (missing(SK)) stop("must supply SK")
    n <- length(SK)
    rval <- .C("wrap_gsw_SP_from_SK",
               SA=as.double(SK), n=as.integer(n), SP=double(n), NAOK=TRUE, package="gsw")$SP
    if (is.matrix(SK))
        dim(rval) <- dim(SK)
    rval
}

#' Calculate Practical Salinity from Reference Salinity
#'
#' @param SR Reference Salinity [ g/kg ]
#' @return Practical Salinity (PSS-78) [ unitless ]
#' @examples 
#' library(testthat)
#' SR <- c(34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600)
#' SP <- gsw_SP_from_SR(SR)
#' expect_equal(SP, c(34.386552667080714, 34.564513505458834, 34.696889296869848,
#'                    34.518231743800094, 34.405762086435850, 34.397799632817147))
#' @seealso The reverse is \code{\link{gsw_SR_from_SP}}; also related are \code{\link{gsw_SP_from_SA}}, \code{\link{gsw_SP_from_SK}} and \code{\link{gsw_SP_from_Sstar}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_SR.html}
gsw_SP_from_SR <- function(SR)
{
    if (missing(SR)) stop("must supply SR")
    n <- length(SR)
    rval <- .C("wrap_gsw_SP_from_SR",
               SA=as.double(SR), n=as.integer(n), SP=double(n), NAOK=TRUE, package="gsw")$SP
    if (is.matrix(SR))
        dim(rval) <- dim(SR)
    rval
}

#' Practical Salinity from Preformed Salinity
#' 
#' @param Sstar Preformed Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
#' @param latitude latitude in decimal degrees north [ -90 ... +90 ]
#' @return Practical Salinity (PSS-78) [ unitless ]
#' @examples 
#' library(testthat)
#' Sstar <- c(34.7115, 34.8912, 35.0247, 34.8436, 34.7291, 34.7197)
#' p <- c(         10,      50,     125,     250,     600,    1000)
#' longitude <- 188
#' latitude <- 4
#' SP <- gsw_SP_from_Sstar(Sstar, p, longitude, latitude)
#' expect_equal(SP, c(34.548646570969929, 34.727538423586189, 34.860549501859502,
#'                    34.681006826476434, 34.568065697992346, 34.560023926979518)) 
#' @seealso
#' \code{\link{gsw_Sstar_from_SP}} does the reverse; \code{\link{gsw_SA_from_Sstar}} is similar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_Sstar.html}
gsw_SP_from_Sstar <- function(Sstar, p, longitude, latitude)
{
    if (missing(Sstar)) stop("must supply Sstar")
    l <- argfix(list(Sstar=Sstar, p=p, longitude=longitude, latitude=latitude))
    if (is.null(l$p)) stop("must supply p")
    if (is.null(l$longitude)) stop("must supply longitude")
    if (is.null(l$latitude)) stop("must supply latitude")
    n <- length(Sstar)
    rval <- .C("wrap_gsw_SP_from_Sstar",
               Sstar=as.double(l$Sstar), p=as.double(l$p), longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(Sstar))
        dim(rval) <- dim(Sstar)
    rval
}

#' Calculate Reference Salinity from Practical Salinity
#'
#' @param SP Practical Salinity (PSS-78) [ unitless ]
#' @return Reference Salinity [ g/kg ]
#' @examples 
#' library(testthat)
#' SP <- c(34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600)
#' SR <- gsw_SR_from_SP(SP)
#' expect_equal(SR, c(34.711611927085727, 34.891255045714303, 35.024882197714305,
#'                    34.844535778285724, 34.731002934857159, 34.722965211428587))
#' @seealso The reverse is \code{\link{gsw_SP_from_SR}}; also related are \code{\link{gsw_SP_from_SA}}, \code{\link{gsw_SP_from_SK}} and \code{\link{gsw_SP_from_Sstar}}.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SR_from_SP.html}
gsw_SR_from_SP <- function(SP)
{
    if (missing(SP)) stop("must supply SP")
    n <- length(SP)
    rval <- .C("wrap_gsw_SR_from_SP",
               SP=as.double(SP), n=as.integer(n), SR=double(n), NAOK=TRUE, package="gsw")$SR
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Convert from Absolute Salinity to Preformed Salinity
#'
#' Calculate Preformed Salinity from Absolute Salinity, pressure,
#' longitude, and latitude.
#'
#' If SA is a matrix and if its dimensions correspond to the
#' lengths of longitude and latitude, then the latter are
#' converted to analogous matrices with \code{\link{expand.grid}}.
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return Preformed Salinity [ g/kg ]
#' @examples
#' library(testthat)
#' SA <-   c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' p <-    c(     10,      50,     125,     250,     600,    1000)
#' lat <-  c(      4,       4,       4,       4,       4,       4)
#' long <- c(    188,     188,     188,     188,     188,     188)
#' Sstar <- gsw_Sstar_from_SA(SA,p,long,lat)
#' expect_equal(Sstar, c(34.711575335926490, 34.891138777337822, 35.024705401162166,
#'                       34.843564118358302, 34.729005527604883, 34.719712883389462))    
#' @seealso \code{\link{gsw_SA_from_Sstar}} does the reverse; \code{\link{gsw_Sstar_from_SP}} is similar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_Sstar_from_SA.html}
gsw_Sstar_from_SA <- function(SA, p, longitude, latitude)
{
    if (missing(SA)) stop("must supply SA")
    if (missing(p)) stop("must supply p")
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    ## check for special case that SP is a matrix defined on lon and lat
    if (is.matrix(SA)) {
        dim <- dim(SA)
        if (length(longitude) == dim[1] && length(latitude) == dim[2]) {
            ll <- expand.grid(longitude=as.vector(longitude), latitude=as.vector(latitude))
            longitude <- ll$longitude
            latitude <- ll$latitude
        }
    }
    l <- argfix(list(SA=SA, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_Sstar_from_SA",
               SA=as.double(l$SA), p=as.double(l$p), longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Convert from Practical Salinity to Preformed Salinity
#'
#' Calculate Preformed Salinity from Practical Salinity, pressure,
#' longitude, and latitude.
#'
#' If SP is a matrix and if its dimensions correspond to the
#' lengths of longitude and latitude, then the latter are
#' converted to analogous matrices with \code{\link{expand.grid}}.
#' 
#' @param SP Practical Salinity (PSS-78) [ unitless ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return Preformed Salinity [ g/kg ]
#' @examples
#' library(testthat)
#' SP <-   c(34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600)
#' p <-    c(     10,      50,     125,     250,     600,    1000)
#' lat <-  c(      4,       4,       4,       4,       4,       4)
#' long <- c(    188,     188,     188,     188,     188,     188)
#' Sstar <- gsw_Sstar_from_SP(SP,p,long,lat)
#' expect_equal(Sstar, c(34.711553680880769, 34.891161395333754, 35.024650265047370,
#'                       34.843593141519356, 34.729033995955525, 34.719675962471783))   
#' @seealso \code{\link{gsw_SP_from_Sstar}} does the reverse; \code{\link{gsw_Sstar_from_SA}} is similar.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_Sstar_from_SP.html}
gsw_Sstar_from_SP <- function(SP, p, longitude, latitude)
{
    if (missing(SP)) stop("must supply SP")
    if (missing(p)) stop("must supply p")
    if (missing(longitude)) stop("must supply longitude")
    if (missing(latitude)) stop("must supply latitude")
    ## check for special case that SP is a matrix defined on lon and lat
    if (is.matrix(SP)) {
        dim <- dim(SP)
        if (length(longitude) == dim[1] && length(latitude) == dim[2]) {
            ll <- expand.grid(longitude=as.vector(longitude), latitude=as.vector(latitude))
            longitude <- ll$longitude
            latitude <- ll$latitude
        }
    }
    l <- argfix(list(SP=SP, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_Sstar_from_SP",
               SP=as.double(l$SP), p=as.double(l$p), longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Freezing temperature
#'
#' This uses the C function named \code{gsw_t_freezing_exact}, because the 
#' C function named \code{gsw_t_freezing} does not produce check values that
#' match the Matlab function called \code{gsw_t_freezing} (see references 
#' for those test values).
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param saturation_fraction saturation fraction of dissolved air in seawater
#' @return in-situ freezing temperature (ITS-90) [ deg C ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' p <- c(      10,      50,     125,     250,     600,    1000)
#' saturation_fraction <- 1
#' tf <- gsw_t_freezing(SA, p, saturation_fraction)
#' expect_equal(tf, c(-1.902730710149803, -1.942908619287183, -2.006861069199743,
#'                    -2.090985086875259, -2.351293130342102, -2.660498762776720)) 
#' @seealso \code{\link{gsw_CT_freezing}} is the analogue for Conservative Temperature.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_t_freezing.html}
#' \url{http://www.teos-10.org/pubs/gsw/v3_04/html/gsw_t_freezing_poly.html}
gsw_t_freezing <- function(SA, p, saturation_fraction=1)
{
    l <- argfix(list(SA=SA, p=p, saturation_fraction=saturation_fraction))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_t_freezing",
               SA=as.double(l$SA), p=as.double(l$p), saturation_fraction=as.double(l$saturation_fraction),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' In situ temperature from Conservative Temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return in-situ temperature (ITS-90) [ deg C ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' t <- gsw_t_from_CT(SA, CT, p)
#' expect_equal(t, c(28.785580227725703, 28.432872246163946, 22.810323087627076,
#'                   10.260010752788906, 6.886286301029376, 4.403624452383043))
#' @seealso \code{\link{gsw_CT_from_t}} does the reverse.
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_t_from_CT.html}
gsw_t_from_CT <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_t_from_CT",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Thermobaric coefficient (75-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return thermobaric coefficient wrt Conservative Temperature [ 1/(K Pa) ]
#' @examples 
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' tb <- gsw_thermobaric(SA, CT, p)
#' expect_equal(tb, 1e-11 * c(0.152618598186650, 0.153662896162852, 0.173429325875738,
#'                            0.232810160208414, 0.251984724005424, 0.266660342289558))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_thermobaric.html}
gsw_thermobaric <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_thermobaric",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Turner angle and density ratio
#'
#' This uses the 48-term density equation. The values of Turner Angle
#' Tu and density ratio Rrho are calculated at mid-point pressures, p_mid.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return list containing Tu [ degrees ], Rsubrho [ unitless ], and p_mid [ dbar ]
#' @examples
#' library(testthat)
#' SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
#' p <-  c(     10,      50,     125,     250,     600,    1000)
#' r <- gsw_Turner_Rsubrho(SA, CT, p)
#' expect_equal(r$Tu, c(-2.063858905281147, 41.758435216784427, 47.606966981687535,
#'                      53.710351151706369, 45.527063858211527)) 
#' expect_equal(r$Rsubrho, 100*c(-0.009304335069039, -0.176564834348709, 0.219627771740757,
#'                               0.065271424662002, 1.087044054679743)) 
#' expect_equal(r$p_mid, 100*c(0.300, 0.875, 1.875, 4.250, 8.000))
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_Turner_Rsubrho.html}
gsw_Turner_Rsubrho <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    r <- .C("wrap_gsw_Turner_Rsubrho",
            SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
            n=n, Tu=double(n-1), Rsubrho=double(n-1), p_mid=double(n-1))
    Tu <- r$Tu
    Rsubrho <- r$Rsubrho
    p_mid <- r$p_mid
    if (is.matrix(SA)) {
        stop("gsw_Turner_Rsubrho() cannot handle matrix SA")
        ## dim(Tu) <- dim(SA)
        ## dim(Rsubrho) <- dim(SA)
        ## dim(p_mid) <- dim(SA)
    }
    list(Tu=Tu, Rsubrho=Rsubrho, p_mid=p_mid)
}

#' Height from pressure (75-term equation)
#' 
#' @param p sea pressure [ dbar ]
#' @param latitude latitude in decimal degrees north [ -90 ... +90 ]
#' @return height [ m ]
#' @examples
#' library(testthat)
#' z <- gsw_z_from_p(c(10, 50, 125, 250, 600,1000), 4)
#' expect_equal(z, 1e2*c(-0.099445834469453, -0.497180897012550, -1.242726219409978,
#'                       -2.484700576548589, -5.958253480356214, -9.920919060719987)) 
#' @seealso
#' This is (almost) the reverse of \code{\link{gsw_p_from_z}}
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_z_from_p.html}
gsw_z_from_p <- function(p, latitude)
{
    l <- argfix(list(p=p, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_z_from_p",
               p=as.double(l$p), lat=as.double(l$latitude),
               n=n, rval=double(n), NAOK=TRUE, package="gsw")$rval
    if (is.matrix(p))
        dim(rval) <- dim(p)
    rval
}

