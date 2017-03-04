# gsw

[![Build Status](https://travis-ci.org/TEOS-10/GSW-R.svg?branch=master)](https://travis-ci.org/TEOS-10/GSW-R)
[![codecov](https://codecov.io/gh/TEOS-10/GSW-R/branch/master/graph/badge.svg)](https://codecov.io/gh/TEOS-10/GSW-R)
![RStudio CRAN mirror downloads](http://cranlogs.r-pkg.org/badges/last-month/gsw)
![RStudio CRAN mirror downloads](http://cranlogs.r-pkg.org/badges/last-week/gsw)
![RStudio CRAN mirror downloads](http://cranlogs.r-pkg.org/badges/last-day/gsw)

gsw is an R package that provides a connection to the thermodynamic equation of
seawater as defined at [teos-10.org](http://www.teos-10.org).  An earlier version
of gsw was referenced to TEOS-10 version 3.03, but the present one is
referenced to version 3.05, which is current as of early 2017. See the
[CRAN](https://cran.r-project.org/package=gsw) website for check statistics,
etc.

All the gsw functions reproduce [teos-10.org](http://www.teos-10.org) test
values to a tolerance of 1.5e-8, the default for numerical comparison in R
working on a 64-bit machine. (The tests are part of the package-building
process, as is usual in R.)

Function names match those used specified in the TEOS-10
[documentation](http://www.teos-10.org/pubs/gsw/html/gsw_contents.html).

