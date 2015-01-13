gsw is an R package that provides a connection to the thermodynamic equation of
seawater version 3.03, downloaded from [teos-10.org](http://www.teos-10.org) on
2014-08-25. This package is still in development, and is on track for a planned
submission to CRAN in January.

All the functions provided here reproduce [teos-10.org](http://www.teos-10.org)
test values to a tolerance of 1.5e-8, the default for numerical comparison in R
working on a 64-bit machine. (The tests are part of the package-building
process, as is usual in R.)

Function names match those used specified in the TEOS-10
[documentation](http://www.teos-10.org/pubs/gsw/html/gsw_contents.html).

The plan is to update this package to match publicised changes to GSW.  Users
should be aware that this may involve changes to function names or arguments,
if TEOS-10 alters such things.

