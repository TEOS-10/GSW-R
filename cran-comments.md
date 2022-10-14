# Submission of 1.0-7

The changes since 1.0-6 are as follows.

- Update to GSW-C as of 2021-07-14, github commit
  `bad2c9e4e154597ce563aaaf3ce09b1c52a2ab46`. This does not change any
  existing functions in GSW-R.
- Add `gsw_SP_salinometer()`.
- Add `gsw_o2sol()`.
- Add `gsw_o2sol_SP_pt()`.


# Tests

## Local Tests

* local macOS (12.5 beta, intel) with R-4.2.0 (2022-04-22)
* win-builder (devel and release)

## Github R-CMD-check Action Tests

* macOS-latest (release)
* windows-latest (release)
* ubuntu-latest (devel)
* ubuntu-latest (release)
* ubuntu-latest (oldrel-1)


