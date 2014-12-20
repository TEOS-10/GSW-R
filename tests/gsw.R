library(gsw)
# Test against with values provided on the TEOS-10 website.
 
CT <- gsw_CT_from_t(34.7118, 28.7856, 10)
stopifnot(all.equal(CT, 28.809919826700281))

g <- gsw_grav(c(-90, -60, -30, 0), 0)
stopifnot(all.equal(g, c(9.832186205884799, 9.819178859991149,
                         9.793249257048750, 9.780327000000000)))

SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
CT <- c(28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236)
p <- c(      10,      50,     125,     250,     600,    1000)
latitude <- 4
N2 <- gsw_Nsquared(SA, CT, p, latitude)
stopifnot(all.equal(N2, 1e-3*c(0.060847042791371, 0.235730438897447,
                               0.216590355253073, 0.012935081031687,
                               0.008430222212653)))

rho <- gsw_rho(34.7118, 28.8099, 10)
stopifnot(all.equal(rho, 1021.8404465661))

SA <- gsw_SA_from_SP(34.5487, 10, 188, 4)
stopifnot(all.equal(SA, 34.711778344814114))

SP <- gsw_SP_from_C(34.5487, 28.7856, 10)
stopifnot(all.equal(SP, 20.009869599086951))

z <- gsw_z_from_p(10, 4)
stopifnot(all.equal(z, -9.9445831334188))
