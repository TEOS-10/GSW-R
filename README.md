# gsw

TravisCI test: [![Build Status](https://travis-ci.org/TEOS-10/GSW-R.svg?branch=master)](https://travis-ci.org/TEOS-10/GSW-R)

CRAN: ![RStudio CRAN mirror downloads](http://cranlogs.r-pkg.org/badges/last-month/gsw)
![RStudio CRAN mirror downloads](http://cranlogs.r-pkg.org/badges/last-week/gsw)
![RStudio CRAN mirror downloads](http://cranlogs.r-pkg.org/badges/last-day/gsw)

gsw is an R package that provides a connection to the thermodynamic equation of
seawater version 3.03, downloaded from [teos-10.org](http://www.teos-10.org) on
2014-08-25. This package was first installed on
[CRAN](http://cran.r-project.org/web/packages/gsw/) in January 2015.

All the gsw functions reproduce [teos-10.org](http://www.teos-10.org) test
values to a tolerance of 1.5e-8, the default for numerical comparison in R
working on a 64-bit machine. (The tests are part of the package-building
process, as is usual in R.)

Function names match those used specified in the TEOS-10
[documentation](http://www.teos-10.org/pubs/gsw/html/gsw_contents.html).

The plan is to update this package to match publicized changes to GSW.  Users
should be aware that this may involve changes to function names or arguments,
if TEOS-10 alters such things.


