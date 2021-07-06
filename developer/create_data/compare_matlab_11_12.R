# Compare data in matlab versions 3.06-11 and 3.06-12

library(R.matlab)
m11 <- readMat("~/src/gsw_matlab_v3_06_11/library/gsw_data_v3_0.mat")
str(m11, 1)

m12 <- readMat("~/src/gsw_matlab_v3_06_12/library/gsw_data_v3_0.mat")
str(m12, 1)

