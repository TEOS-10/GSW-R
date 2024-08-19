# Submission of 1.2-0

This update was required by 2024-09-20, in order to accommodate the
CRAN-requested change (on that date) relating to the need to replace the
Calloc() and Free() calls with R_Calloc() and R_Free() calls. I wish to thank
the CRAN team for sending that email, since I do not check the CRAN-results
pages as frequently as I ought to, for this mainly-stable package.

The only other change from the previous CRAN version is the addition of the
function `gsw_infunnel()`.  This is mainly for developers' use and, in any
case, its provision does not affect the actions of the rest of the package, so
it should not affect users adversely.

# Tests

## Local Tests

I saw no problems with R 4.4.1 on macOS Beta 15.0 Beta (24A5298h) revealed no
ERRORs, no WARNINGs, or NOTEs (apart from the usual NOTE about the author
name).

## Remote tests

* I saw no problems on win-builder (devel and release).
* I saw no problems on Github R-CMD-check action tests.
