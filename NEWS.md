# gsw 1.2-0

* add `gsw_infunnel()`
* switch from `Calloc` to `R_Calloc` etc, so package builds with "strict"
  header checking (issue #63).

# gsw 1.1-1

* Remove a journal link that fails because of an invalid certificate.
* Fix a C function declaration that is invalid for some compilers.

# gsw 1.1-0

This update is based on C code from <https://github.com/TEOS-10/GSW-C> (release
v3.06-16-0, commit 657216dd4f5ea079b5f0e021a4163e2d26893371) and the .mat data
file from <https://github.com/TEOS-10/GSW-Matlab> (commit
38c9635d6fd93e74c2648e4ee23cec49c1f58530). It was necessary to make some
changes to the within-documentation test suite, as listed below; MRD stands for
Mean Relative Difference from the values in the previous CRAN release of `gsw`,
version 1.0-6 released 2021-07-07.

* `gsw_rho_second_derivatives()$rho_SA_p`: MRD 0.4
* `gsw_rho_second_derivatives_wrt_enthalpy()$rho_SA_SA`: MRD 1e-3
* `gsw_specvol_second_derivatives()$specvol_SA_p`: MRD 0.2
* `gsw_specvol_second_derivatives_wrt_enthalpy()$specvol_SA_SA`: MRD 3e-4
* `gsw_thermobaric()$rb`: MRD 0.04

# gsw 1.0-7

* Update to GSW-C as of 2021-07-14, github commit
  `bad2c9e4e154597ce563aaaf3ce09b1c52a2ab46`. This does not change any
  existing functions in GSW-R.
* Add gsw_SP_salinometer().
* Add gsw_o2sol().
* Add gsw_o2sol_SP_pt().

# gsw 1.0-6

* Update to GSW-C as of 2021-07-06, github commit
  `9c10670e89fce906da2cebce3399d73c054e769e`.
* Remove dependency on the 'testthat' package.
* gsw_z_from_p() and gsw_p_from_z() gain parameters geo_strf_dyn_height and
  sea_surface_geopotential.
* Relax gsw_geo_strf_dyn_height() test, since values are tied to a new Matlab
  interpolation scheme that is not incorporated in GSW-C.

# gsw 1.0-5

* Update to GSW-C 3.05-4.
* Make internal sort routines be machine-independent.

# gsw 1.0-4

* Update to GSW-C 3.05-3.
* Handle NaN values better.
* Add dozens of new GSW functions.

# gsw 1.0-3

* First version, using GSW version 3.03.
