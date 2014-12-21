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
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_alpha.html}
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

#' Conductivity from Practical Salinity
#' 
#' @param SP Practical Salinity (PSS-78) [ unitless ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @examples 
#' gsw_C_from_SP(34.5487, 28.7856, 10) # 56.412599581571186
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_C_from_SP.html}
gsw_C_from_SP <- function(SP, t, p)
{
    l <- argfix(list(SP=SP, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_C_from_SP",
               SP=as.double(l$SP), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n))$rval
    if (is.matrix(SP))
        dim(rval) <- dim(SP)
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
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_beta.html}
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

#' Gravitational acceleration
#' 
#' @param latitude latitude in decimal degress north [ -90 ... +90 ]
#' @param p sea pressure [ dbar ]
#' @return gravitational acceleration [ m/s^2 ]
#' @examples
#' gsw_grav(c(-90, -60), 0) # 9.832186205884799, 9.819178859991149
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_grav.html}
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

#' Calculate Brunt Vaisala Frequency squared
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @param latitude latitude in decimal degrees [ -90 to 90 ]
#' @return a list containing N2 [ s^(-2) ] and mid-point pressure p_mid [ dbar ]
#' @examples 
#' SA <- c(34.7118, 34.8915)
#' CT <- c(28.8099, 28.4392)
#' p <- c(      10,      50)
#' latitude <- 4
#' gsw_Nsquared(SA, CT, p, latitude)$N2 # 6.0847042791371e-5
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html}
gsw_Nsquared <- function(SA, CT, p, latitude=0)
{
    l <- argfix(list(SA=SA, CT=CT, p=p, latitude=latitude))
    n <- length(l[[1]])
    r <- .C("wrap_gsw_Nsquared",
            SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p), latitude=as.double(l$latitude),
            n=n, n2=double(n-1), p_mid=double(n-1))
    N2 <- r$n2
    p_mid <- r$p_mid
    ## How to handle the reduction in length for a matrix??
    if (is.matrix(SA)) {
        stop("gsw_Nsquared() cannot handle matix SA")
    }
    list(N2=N2, p_mid=p_mid)
}

#' potential density
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @param p_ref reference pressure [ dbar ]
#' @return potential density [ kg/m^3 ]
#' @examples
#' gsw_pot_rho_t_exact(34.7118, 28.7856, 10, 0) # 1021.798145811089
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_pot_rho_t_exact.html}
gsw_pot_rho_t_exact <- function(SA, t, p, p_ref)
{
    l <- argfix(list(SA=SA, t=t, p=p, p_ref=p_ref))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_pot_rho_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p), pref=as.double(l$p_ref),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
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
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_rho.html}
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

#' in-situ density
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return in-situ density [ kg/m^3 ]
#' @examples
#' gsw_rho_t_exact(34.7118, 28.7856, 10) # 1021.840173185531
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_rho_t_exact.html}
gsw_rho_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_rho_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
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

#' potential density anomaly referenced to 0 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 0 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly [ kg/m^3 ]
#' @examples
#' gsw_sigma0(34.7118, 28.8099) # 21.798411276610750
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma0.html}
gsw_sigma0 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma0",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' potential density anomaly referenced to 1000 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 1000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly [ kg/m^3 ]
#' @examples
#' gsw_sigma1(34.7118, 28.8099) # 25.955891533636986
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma1.html}
gsw_sigma1 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma1",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' potential density anomaly referenced to 2000 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 2000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly [ kg/m^3 ]
#' @examples
#' gsw_sigma2(34.7118, 28.8099) # 30.022796416066058
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma2.html}
gsw_sigma2 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma2",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' potential density anomaly referenced to 3000 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 3000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly with reference pressure 3000 dbar [ kg/m^3 ]
#' @examples
#' gsw_sigma3(34.7118, 28.8099) # 34.002600253012133
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma3.html}
gsw_sigma3 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma3",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' potential density anomaly referenced to 4000 dbar
#'
#' This uses the 48-term density equation, and returns
#' potential density referenced to a pressure of 4000 dbar,
#' minus 1000 kg/m^3.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @return potential density anomaly with reference pressure 4000 dbar [ kg/m^3 ]
#' @examples
#' gsw_sigma3(34.7118, 28.8099) # 37.898467323406976
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sigma4.html}
gsw_sigma4 <- function(SA, CT)
{
    l <- argfix(list(SA=SA, CT=CT))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sigma4",
               SA=as.double(l$SA), CT=as.double(l$CT),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' sound speed with 48-term density
#'
#' This uses the 48-term density equation.
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param CT Conservative Temperature [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return sound speed [ m/s ]
#' @examples
#' gsw_sound_speed(34.7118, 28.7856, 10) # 1542.420534932182
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sound_speed.html}
gsw_sound_speed<- function(SA, CT, p)
{
    l <- argfix(list(SA=SA, CT=CT, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sound_speed",
               SA=as.double(l$SA), CT=as.double(l$CT), p=as.double(l$p),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
    rval
}

#' sound speed
#'
#' @param SA Absolute Salinity [ g/kg ]
#' @param t in-situ temperature (ITS-90) [ deg C ]
#' @param p sea pressure [ dbar ]
#' @return sound speed [ m/s ]
#' @examples
#' gsw_sound_speed_t_exact(34.7118, 28.7856, 10) # 1542.420534932182
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_sound_speed_t_exact.html}
gsw_sound_speed_t_exact <- function(SA, t, p)
{
    l <- argfix(list(SA=SA, t=t, p=p))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_sound_speed_t_exact",
               SA=as.double(l$SA), t=as.double(l$t), p=as.double(l$p),
               n=n, rval=double(n))$rval
    if (is.matrix(SA))
        dim(rval) <- dim(SA)
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

#' Freezing temperature
#' 
#' @param SA Absolute Salinity [ g/kg ]
#' @param p sea pressure [ dbar ]
#' @param saturation_fraction saturation fraction of dissolved air in seawater
#' @return in-situ freezing temperature (ITS-90) [ deg C ]
#' @examples 
#' gsw_t_freezing(34.7118, 10) # -1.902730710149803
#' @references
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_t_freezing.html}
gsw_t_freezing <- function(SA, p, saturation_fraction=1)
{
    l <- argfix(list(SA=SA, p=p, saturation_fraction=saturation_fraction))
    n <- length(l[[1]])
    rval <- .C("wrap_gsw_t_freezing",
               SA=as.double(l$SA), p=as.double(l$p), saturation_fraction=as.double(l$saturation_fraction),
               n=n, rval=double(n))$rval
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
#' @return a list containing Tu, Rrho, and p_mid
#' @examples
#' SA = c(34.7118, 34.8915)
#' CT = c(28.8099, 28.4392)
#' p =  c(     10,      50)
#' r <- gsw_Turner_Rsubrho(SA, CT, p) # -2.064830032393999, -0.9304018848608, 30
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
        stop("gsw_Turner_Rsubrho() cannot handle matix SA")
        ## dim(Tu) <- dim(SA)
        ## dim(Rsubrho) <- dim(SA)
        ## dim(p_mid) <- dim(SA)
    }
    list(Tu=Tu, Rsubrho=Rsubrho, p_mid=p_mid)
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
#' \url{http://www.teos-10.org/pubs/gsw/html/gsw_z_from_p.html}
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

