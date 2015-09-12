library(gsw)
## Test against with values provided on the TEOS-10 website, e.g.
##   http://www.teos-10.org/pubs/gsw/html/gsw_adiabatic_lapse_rate_from_t.html
## for the first function.

## gsw_adiabatic_lapse_rate_from_CT()
SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
CT <- c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
p <- c(      10,      50,     125,     250,     600,    1000)
lr <- gsw_adiabatic_lapse_rate_from_CT(SA, CT, p)
expect_equal(lr, 1e-7*c(0.240199646230069, 0.238457486976761, 0.203635157319712,
                        0.119829566859790, 0.100052760967308, 0.087773070307283))

## gsw_alpha()
SA <- c(34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324)
CT <- c(28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036)
p <-  c(     10,      50,     125,     250,     600,    1000)
alpha <- gsw_alpha(SA,CT,p)
expect_equal(alpha, 1e-3 * c( 0.324464211877393, 0.322610094680523, 0.281335030247435,
                             0.173529986885424, 0.146898108553385, 0.130265123640082))
