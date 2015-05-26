How to add new functions
========================

**Preamble.** There are a great many functions in the GSW library, and adding
each one to this package takes a fair bit of time. The steps for adding
functions are given in detail here, in hopes that some users will contribute to
the coding process. Please send contributions to the developers with ``git
pull`` requests.

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

1. Write a C-language wrapper (in ``src/wrappers.c``) for the new GSW function
   you wish to add.  The function name should match that in the [matlab
library](http://www.teos-10.org/pubs/gsw/html/gsw_contents.html), including the
case notation, but note that the GSW C functions that you will be calling all
all named in *lower* case. It should not be difficult for someone with moderate
C skill to write a wrapper; just examine the existing code and be aware of
whether the underlying GSW function returns ``void`` or a value, and mimic a
wrapper of similar type. You may test your C code for compilation errors by
executing

        R CMD SHLIB wrappers.c
in the ``src`` directory.

2. Once a C wrapper exists, create an R function that hooks up with it. Again,
   use a function name that exactly matches the in the [matlab
library](http://www.teos-10.org/pubs/gsw/html/gsw_contents.html), particularly
with regard to case. Any reasonable R programmer should be able to write a
wrapper easily, by mimicking an existing function. Be sure to copy and modify
the reshaping code; without this, many circumstances will yield C problems
(e.g. going past the end of a vector) that will cause R to crash. Be sure to
document the function within ``R/gsw.R`` using the Roxygen notation, mimicking
the other documentation. And be sure to use the capitalization used on the
TEOS-10 website for matlab code, because that is what users will probably
expect. Place the function in its correct alphabetical position.

3. Once the function is added, execute

        Rscript -e "roxygen2::roxygenise()"
in the command line, to rebuild the documentation.

4. Add an entry for the function in ``../NAMESPACE``.

5. Add a test to ``../tests/gsw.R``, using the test values in the GSW
   documentation provided on the TEOS-10 website.

6. Build the package to see that code matches docs, etc. This can be done in
   Rstudio, or in the commandline with e.g.

        R CMD BUILD GSW-R
        R CMD CHECK gsw_0.1-1.tar.gz
   where of course the version number may need to be change.

