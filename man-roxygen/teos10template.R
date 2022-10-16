#' @section Implementation Note:
#'
#' This R function uses a wrapper to a C function contained within the GSW-C
#' system as updated 2022-10-11 at \url{https://github.com/TEOS-10/GSW-C} with
#' git commit `657216dd4f5ea079b5f0e021a4163e2d26893371`.
#'
#' The C function uses data from the \code{library/gsw_data_v3_0.mat}
#' file provided in the GSW-Matlab source code, version 3.06-11.
#' Unfortunately, this version of the mat file is no longer displayed on the
#' TEOS-10.org website.  Therefore, in the interests of making GSW-R be
#' self-contained, a copy was downloaded from
#' \url{http://www.teos-10.org/software/gsw_matlab_v3_06_11.zip} on 2022-05-25,
#' the .mat file was stored in the developer/create_data directory of
#' \url{https://github.com/TEOS-10/GSW-R}, and then the dataset used in GSW-R
#' was created based on that .mat file.
#'
#' Please consult \url{http://www.teos-10.org} to learn more about the various
#' TEOS-10 software systems.

