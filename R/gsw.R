## Please see ../README_developer.md for the scheme used in adding functions
## here. Generally the functions will be added to Part 4.


## PART 1: document the package

#' R implementation of the Thermodynamic Equation Of Seawater - 2010 (TEOS-10)
#'
#' Provides an R interface to the TEOS-10 / GSW (Gibbs Sea Water) library.
#' The functions are designed to match the Matlab implementation, so the
#' function names match and the documentation for each function in the 
#' present package contains a link to the official TEOS-10 documentation
#' of the paired Matlab function.
#'
#' As of late 2014, the package is still in an early stage of development,
#' with only a handful of (important) functions working.
#'
#' @docType package
#' @name gsw
NULL


## PART 2: utility functions

#' Reshape list elements to match the shape of the first element.
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

#' thermal expansion coefficient with respect to Conservative Temperature. (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return thermal expansion coefficient with respect to Conservative Temperature [ 1/K ]
#' @examples
#' SA = c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT = c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p =  c(     10,      50,     125,     250,     600,    1000)
#' alpha <- gsw_alpha(SA,CT,p)
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_SP.html}
gsw_alpha <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_alpha",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' saline contraction coefficient at constant Conservative Temperature. (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return saline contraction coefficient at constant Conservative Temperature [ kg/g ]
#' @examples
#' SA = c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
#' CT = c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
#' p =  c(     10,      50,     125,     250,     600,    1000)
#' beta <- gsw_beta(SA,CT,p)
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_SP.html}
gsw_beta <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_beta",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Convert from temperature to conservative temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @examples 
#' gsw_CT_from_t(34.7118, 28.7856, 10) # 28.809919826700281
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_CT_from_t.html}
gsw_CT_from_t <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_CT_from_t",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n))$rval
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
#' @return Brunt-Vaisala Frequency squared [ s^(-2) ]
#' @examples 
#' SA <- c(34.7118, 34.8915)
#' CT <- c(28.8099, 28.4392)
#' p <- c(      10,      50)
#' latitude <- 4
#' gsw_Nsquared(SA, CT, p, latitude) # 6.0847042791371e-5
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html}
gsw_Nsquared <- function(SA, CT, p, latitude=0)
{
    l <- argfix(list(SA=SA, CT=CT, p=p, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_Nsquared",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p), latitude=as.double(l$latitude),
               n=n, n2=double(n), p_mid=double(n))$n2
    ## How to handle the reduction in length for a matrix??
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    else
        rval <- head(rval, -1)
    rval
}

#' Gravitational acceleration
#' 
#' @param latitude latitude in decimal degress north [ -90 ... +90 ]
#' @param p sea pressure [ dbar ]
#' @return gravitational acceleration [ m/s^2 ]
#' @examples
#' gsw_grav(c(-90, -60, -30, 0), 0) # 9.832186205884799, 9.819178859991149, 9.793249257048750, 9.780327000000000
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_SP.html}
gsw_grav <- function(latitude, p)
{
    l <- argfix(list(latitude=latitude, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_grav",
               latitude=as.double(l$latitude), p=as.double(l$p),
               n=n, rval=double(n))$rval
    if (is.matrix(latitude))
        dim(rval) <- dim(latitude)
    rval
}

#' in-situ density (48-term equation)
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return in-situ density [ kg/m^3 ]
#' @examples
#' gsw_rho(34.7118, 28.8099, 10) # 1021.8404465661
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_SP.html}
gsw_rho <- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_rho",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' Convert from practical salinity to absolute salinity
#' 
#' @param SP Practical Salinity (PSS-78) [ unitless ]
#' @param p sea pressure [ dbar ]
#' @param longitude longitude in decimal degrees [ 0 to 360 or -180 to 180]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @examples
#' gsw_SA_from_SP(34.5487, 10, 188, 4) # 34.711778344814114 
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_SP.html}
gsw_SA_from_SP <- function(SP, p, longitude, latitude)
{
    l <- argfix(list(SP=SP, p=p, longitude=longitude, latitude=latitude))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SA_from_SP",
               SA=as.double(l$SP), p=as.double(l$p), longitude=as.double(l$longitude), latitude=as.double(l$latitude),
               n=n, rval=double(n))$rval
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
    rval
}

#' Convert from conductivity to practical salinity
#' 
#' @param C conductivity [ mS/cm ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return Practical salinity.
#' @examples 
#' gsw_SP_from_C(34.5487, 28.7856, 10) # 20.009869599086951
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_C.html}
gsw_SP_from_C <- function(C, t, p)
{
    l <- argfix(list(C=C, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_SP_from_C",
               C=as.double(l$C), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n))$rval
    if (is.matrix(C))
        dim(rval) <- dim(C)
    rval
}

#' height from pressure (48-term equation)
#' 
#' @param p sea pressure [ dbar ]
#' @param lat latitude in decimal degrees north [ -90 ... +90 ]
#' 
#' @return height [ m ]
#' @examples
#' gsw_z_from_p(10, 4) # -9.9445831334188
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_SP.html}
gsw_z_from_p<- function(p, lat)
{
    l <- argfix(list(p=p, lat=lat))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_z_from_p",
               p=as.double(l$p), lat=as.double(l$lat),
               n=n, rval=double(n))$rval
    if (is.matrix(p))
        dim(rval) <- dim(p)
    rval
}

