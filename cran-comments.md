# Submission of 1.2-0

This update was required by 2024-09-20, in order to accommodate the change from
the Calloc() and Free() macros to R_Calloc() and R_Free(). The only other
change from the previous CRAN version is the addition of the function
`gsw_infunnel()`, and this addition does not affect any previous behaviours of
the package.

# Tests

## Local Tests

R 4.4.1 on macOS Beta 15.0 Beta (24A5298h) revealed no ERRORs, no WARNINGs,
or NOTEs. One NOTE was about the author name and the other about the package
size (a known issue; see above).  These tests included those in the CRAN test
suite, along with other tests of datasets that are not provided with the
package.

* No problems on macOS (12.5 beta, intel) with R-4.2.1 (2022-06-23).

## Remote tests

* No problems on win-builder (devel and release).
* No problems on Github R-CMD-check action tests.
