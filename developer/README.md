How to update for a new version of GSW
======================================

The C code that provides the foundation for this is based on the [C
version](https://github.com/TEOS-10/GSW-C), which I think is mainly based on
the [Fortran version](https://github.com/TEOS-10/GSW-Fortran).

Updating the gsw package requires the following steps.


1. Use `git` to update local copies of both the C and Fortran versions, and
   store them somewhere convenient (in the illustrations provided here, they
are in `~/git`, )
2. Copy the C version as follows

    ```
    cd src # from the main 'gsw' directory
    cp  ~/git/GSW-C/gsw_internal_const.h .
    cp  ~/git/GSW-C/gsw_oceanographic_toolbox.c gsw_oceanographic_toolbox.c
    cp  ~/git/GSW-C/gswteos-10.h gswteos-10.h 
    ```
3. Look in `man/saar.Rd` to find the md5 value for the netcdf file that contains the SA atlas data. Then check that value against the value of the file named in `developer/create_data/create_data.R`. If the values differ, that means that the SA atlas file has been updated, and you will need to perform the next step. If they are the same, skip the next step.

4. Depending on the results of the previous step, either skip this step or continue with it by entering `create_data` and doing

    ```
R --no-save < create_data.R
```
Look at the output, and if everything worked (if you don't know how to tell, stop this exercise now) and do.

    ```
cp saar.rda ../data
```
to update the data file.

5. Do a build-check, e.g. in Rstudio or otherwise. If the build fails because of a missing compiler, stop now, because these instructions are meant for developers, all of whom are certain to have compilers installed.

5. Try a test check. This should pass with no errors, no warnings, and no notes. Otherwise, please open an issue.


How to update the vignette
==========================

It is easy to change the text, but if figures are to be changed, things are complicated.

1. Edit `vignettes/gsw.Rmd`, *not* `inst/doc/gsw.Rmd`. (The latter is created from the former by the building process, so if you edit it, your changes will get clobbered.)

2. If your edits do not involve diagrams, run `devtools::build_vignettes()` to update the vignette. This creates some new files and makes copies as appropriate. Do not do alter the resultant files because, again, your changes will get clobbered upon another build.

3. If your edits involve diagrams, a trick must be done. This is because the diagrams are made with the `oce` package, but `oce` depends on `gsw`, and thus `devtools::build_vignettes()` would get into a circular condition if it tried to update the diagrams. For this reason, diagram-producing code within the vignette is not actually run. (If you cannot see that fact from examining `vignettes/gsw.Rmd`, you should probably not be trying to make changes and should contact the authors instead with suggested revisions. The procedure for making new diagrams is quite simple. Open an R session in the `vignettes` directory, and run `knitr::purl("gsw.Rmd")`. That will create a file named `gsw.R`.  Run that through R, and you will get some new graphics files. Then, return to R and run `devtools::build_vignettes()` and you'll find that the vignette has been updated with new diagrams. (Note that this extra `purl` step is not commonly required for R packages.)


How to add new functions
========================

**Preamble.** There are many functions in the GSW library, and adding each one
to this package takes a fair bit of time. The steps for adding functions are
given in detail here, in hopes that some users will contribute to the coding
process. Please send contributions to the developers with ``git pull``
requests.

This package performs its work through linkage with TEOS-10 functions written
in the C language. The R function ``.C()`` is used for this linkage.  Because
``.C()`` cannot work with C functions that return a value, it is necessary to
write C-language wrappers for such functions. These wrappers are stored in
``src/wrappers.c``. Some data reshaping is commonly required before calling
these C wrappers, and this is done in with R functions in the file ``R/gsw.R``.

The detailed steps in adding these C/R function pairs are given below. Each
step is important, and if you stumble, you might find it easier to simply open
an issue on the website, explaining to others why the function you need should
be put higher on the priority list.

You will need reasonable R and C skills to do this work. For example, `R/gsw.R`
contains a function called `argfix`. If you cannot understand what that does
within a couple of seconds, you should not be trying to add a new function to
gsw.  Another test of your skill is whether can understand what is being done
in `src/wrappers.c` within a couple of minutes of reading. Unless you are quite
confident in your ability to understand how gsw is constructed, and willing to
following the coding conventions, your best plan would be to open an issue
report and let the developers do the work.

With those provisos, here is a summary of how to add new functions to gsw.

1. Write a C-language wrapper (in ``src/wrappers.c``) for the new GSW function
   you wish to add.  The function name should match that in the [matlab
library](http://www.teos-10.org/pubs/gsw/html/gsw_contents.html), including the
case notation, but note that the GSW C functions that you will be calling all
all named in *lower* case. It should not be difficult for someone with moderate
C skill to write a wrapper; just examine the existing code and be aware of
whether the underlying GSW function returns ``void`` or a value, and mimic a
wrapper of similar type.

2. Once a C wrapper exists, create an R function that hooks up with it. Again,
   use a function name that matches one in the [matlab
library](http://www.teos-10.org/pubs/gsw/html/gsw_contents.html), so users will
be able to find it easily. Place your function alphabetically within `R/gsw.R`
(after the preamble functions in that file.) Follow the coding scheme used by
existing functions, with respect to indentation, variable names, etc. Be sure
to understand what `argfix` does ... if not 

. You must include Roxygen documentation. Use templates for
any arguments that are found in the rest of the code -- to not make up new
`@param` text if you can avoid it, because we want uniformity throughout the
package. Include test values from the teos-10 HTML documentation, and to link
to that webpage in the Reference section.

3. Once the function is added, execute

        Rscript -e "roxygen2::roxygenise()"
in the command line, to rebuild the documentation.

4. Add an entry for the function in ``../NAMESPACE``.

5. Add a test to ``../tests/gsw.R``, using the test values in the GSW
   documentation provided on the TEOS-10 website.

6. Build the package to see that code matches docs, etc. This can be done in
   Rstudio, or in the commandline. The developers use a Makefile that contains as follows

    GSWVSN=$(shell awk '/Version/{print($$2)}' GSW-R/DESCRIPTION)
    gsw: force
            cd GSW-R ; echo "devtools::document(roclets=c('rd', 'collate', 'vignette'))" | R --no-save
            PKG_CFLAGS=--pedantic R CMD BUILD --compact-vignettes="gs+qpdf" GSW-R
            R CMD CHECK --as-cran gsw_${GSWVSN}.tar.gz
            R CMD INSTALL gsw_${GSWVSN}.tar.gz

