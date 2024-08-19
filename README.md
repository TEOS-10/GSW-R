# gsw

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/gsw)](https://cran.r-project.org/package=gsw)
[![R-CMD-check](https://github.com/TEOS-10/GSW-R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TEOS-10/GSW-R/actions/workflows/R-CMD-check.yaml)
[![R-hub](https://github.com/TEOS-10/GSW-R/actions/workflows/rhub.yaml/badge.svg)](https://github.com/TEOS-10/GSW-R/actions/workflows/rhub.yaml)
![RStudio CRAN mirror downloads](https://cranlogs.r-pkg.org/badges/last-month/gsw)
![RStudio CRAN mirror downloads](https://cranlogs.r-pkg.org/badges/last-week/gsw)
![RStudio CRAN mirror downloads](https://cranlogs.r-pkg.org/badges/last-day/gsw)
![RStudio CRAN mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/gsw)

gsw is an R package that provides a connection to software relating to TEOS,
the Thermodynamic Equation Of Seawater (see <http://www.teos-10.org>). This
connection involves R wrappers to C functions within the GSW-C library
(<https://github.com/TEOS-10/GSW-C>, release v3.06-16-0, commit
657216dd4f5ea079b5f0e021a4163e2d26893371), along with a data file within the
GSW-Matlab library
(<https://github.com/TEOS-10/GSW-Matlab/blob/master/Toolbox/library/gsw_data_v3_0.mat>,
commit 38c9635d6fd93e74c2648e4ee23cec49c1f58530).

The foundational algorithms upon which both GSW-C and GSW-R rest were devised
by the Scientific Committee on Oceanic Research / International Association for
the Physical Sciences of the Oceans, Working Group 127.  These algorithms were
recommended as a replacement for a previous system (called UNESCO-80) at the
twenty-fifth assembly in the Intergovernmental Oceanographic Commission in
2009.

The gsw functions reproduce test values in the GSW-Matlab documentation to a
tolerance of 1.5e-8, the default for numerical comparison in R working on a
64-bit machine.  This offers some assurance that the coding process has not
introduced errors.

The naming convention in gsw is patterned on the system described at
<http://www.teos-10.org/pubs/gsw/html/gsw_contents.html>, as a way to reduce
the effort users might face in switching to R.
