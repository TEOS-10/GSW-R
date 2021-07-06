In GSW-R 1.0-5, the `data/saar.rda` file was built up from a NetCDF file
provided by GSW-Fortran.  However, I learned from Eric Firing that the GSW-C
code was built instead on the matlab file provided in the source for GSW
3.06-11, and so that is what is used in this version of GSW-R.

Unfortunately, the TEOS-10 website no longer links to 3.06-11, but it *is*
available at http://www.teos-10.org/software/gsw_matlab_v3_06_11.zip as of
2021-jul-6.  To get around worries that this link might fail in the future, I
have downloaded this source zipfile, expanded it, and copied the relevant
matlab file to this directory, using

    cp /Users/kelley/src/gsw_matlab_v3_06_11/library/gsw_data_v3_0.mat .

and this is what I'll use to build up the .rda file.


