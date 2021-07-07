1.0-6
- Update to GSW-C commit 9c10670e89fce906da2cebce3399d73c054e769e (2021-07-06).
- Remove dependence on the 'testthat' package.
- gsw_z_from_p() and gsw_p_from_z() gain parameters geo_strf_dyn_height and
  sea_surface_geopotential.
- Relax gsw_geo_strf_dyn_height() test, since values are tied to the Matlab
  interpolation scheme, not incorporated in GSW-C.

1.0-5
- Update to GSW-C 3.05-4.
- Make internal sort routines be machine-independent.

1.0-4
- Update to GSW-C 3.05-3.
- Handle NaN values better.
- Add dozens of new GSW functions.

1.0-3
- First version, using GSW version 3.03.

