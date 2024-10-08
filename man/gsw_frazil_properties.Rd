% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsw.R
\name{gsw_frazil_properties}
\alias{gsw_frazil_properties}
\title{Properties of Frazil ice}
\usage{
gsw_frazil_properties(SA_bulk, h_bulk, p)
}
\arguments{
\item{SA_bulk}{Absolute Salinity of a combination of seawater and ice [ g/kg ]}

\item{h_bulk}{enthalpy of a mixture of seawater and ice [ J/kg ]}

\item{p}{sea pressure [dbar], i.e. absolute pressure [dbar] minus 10.1325 dbar}
}
\value{
a list containing \code{SA_final}, \code{h_final} and \code{w_Ih_final}.
}
\description{
Calculation of Absolute Salinity, Conservative Temperature, and ice mass fraction
based on bulk Absolute Salinity, bulk enthalpy, and pressure
}
\section{Implementation Note}{


This R function uses a wrapper to a C function contained within the GSW-C
system as updated 2022-10-11 at \url{https://github.com/TEOS-10/GSW-C} with
git commit `657216dd4f5ea079b5f0e021a4163e2d26893371`.

The C function uses data from the \code{library/gsw_data_v3_0.mat}
file provided in the GSW-Matlab source code, version 3.06-11.
Unfortunately, this version of the mat file is no longer displayed on the
TEOS-10.org website.  Therefore, in the interests of making GSW-R be
self-contained, a copy was downloaded from
\url{http://www.teos-10.org/software/gsw_matlab_v3_06_11.zip} on 2022-05-25,
the .mat file was stored in the developer/create_data directory of
\url{https://github.com/TEOS-10/GSW-R}, and then the dataset used in GSW-R
was created based on that .mat file.

Please consult \url{http://www.teos-10.org} to learn more about the various
TEOS-10 software systems.
}

\examples{
SA_bulk <- c(  34.7118,   34.8915,   35.0256,   34.8472,   34.7366,   34.7324)
h_bulk <- c( -4.5544e4, -4.6033e4, -4.5830e4, -4.5589e4, -4.4948e4, -4.4027e4)
p <- c(             10,        50,       125,       250,       600,      1000)
r <- gsw_frazil_properties(SA_bulk, h_bulk, p)
stopifnot(all.equal(r$SA_final, c(39.111030663000442, 39.407625769681573, 39.595789974885108,
                                  39.481230045372889, 39.591177095552503, 39.826467709177123)))
stopifnot(all.equal(r$CT_final, c(-2.156311126114311, -2.204672298963783, -2.273689262333450,
                                  -2.363714136353600, -2.644541000680772, -2.977651291726651)))
stopifnot(all.equal(r$w_Ih_final, c(0.112480560814322, 0.114600300867556, 0.115421108602301,
                                    0.117372990660305, 0.122617649983886, 0.127906590822347)))
}
\references{
\url{http://www.teos-10.org/pubs/gsw/html/gsw_frazil_properties.html}
}
