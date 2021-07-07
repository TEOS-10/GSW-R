pkgname <- "seacarb"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('seacarb')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("K0")
### * K0

flush(stderr()); flush(stdout())

### Name: K0
### Title: Henry's constant mol/(kg/atm)
### Aliases: K0
### Keywords: utilities

### ** Examples

  K0(S=35,T=25,P=0)



cleanEx()
nameEx("K1")
### * K1

flush(stderr()); flush(stdout())

### Name: K1
### Title: First dissociation constant of carbonic acid (mol/kg)
### Aliases: K1
### Keywords: utilities

### ** Examples

  K1(S=35,T=25,P=0,k1k2="l",pHscale="T")



cleanEx()
nameEx("K1p")
### * K1p

flush(stderr()); flush(stdout())

### Name: K1p
### Title: First dissociation constant of phosphoric acid (mol/kg)
### Aliases: K1p
### Keywords: utilities

### ** Examples

  K1p(35,25,0)



cleanEx()
nameEx("K2")
### * K2

flush(stderr()); flush(stdout())

### Name: K2
### Title: Second dissociation constant of carbonic acid (mol/kg)
### Aliases: K2
### Keywords: utilities

### ** Examples

  K2(35,25,0)



cleanEx()
nameEx("K2p")
### * K2p

flush(stderr()); flush(stdout())

### Name: K2p
### Title: Second dissociation constant of phosphoric acid (mol/kg)
### Aliases: K2p
### Keywords: utilities

### ** Examples

  K2p(35,25,0)



cleanEx()
nameEx("K2si")
### * K2si

flush(stderr()); flush(stdout())

### Name: K2si
### Title: Second dissociation constant of Si(OH)4
### Aliases: K2si
### Keywords: utilities

### ** Examples

  K2si(S=35, T=25, P=0, pHscale="T")



cleanEx()
nameEx("K3p")
### * K3p

flush(stderr()); flush(stdout())

### Name: K3p
### Title: Third dissociation constant of phosphoric acid (mol/kg)
### Aliases: K3p
### Keywords: utilities

### ** Examples

  K3p(35,25,0)



cleanEx()
nameEx("Kb")
### * Kb

flush(stderr()); flush(stdout())

### Name: Kb
### Title: Dissociation constant of boric acid (mol/kg)
### Aliases: Kb
### Keywords: utilities

### ** Examples

  Kb(S=35,T=25,P=0,pHscale="T")



cleanEx()
nameEx("Kf")
### * Kf

flush(stderr()); flush(stdout())

### Name: Kf
### Title: Equilibrium constant of hydrogen fluoride (mol/kg)
### Aliases: Kf
### Keywords: utilities

### ** Examples

  Kf(S=35,T=25,P=0,kf="pf",pHscale="T")



cleanEx()
nameEx("Khs")
### * Khs

flush(stderr()); flush(stdout())

### Name: Khs
### Title: Dissociation constant of hydrogen sulfide (mol/kg)
### Aliases: Khs
### Keywords: utilities

### ** Examples

  Khs(S=35,T=25,P=0, pHscale="T")



cleanEx()
nameEx("Kn")
### * Kn

flush(stderr()); flush(stdout())

### Name: Kn
### Title: Dissociation constant of ammonium (mol/kg)
### Aliases: Kn
### Keywords: utilities

### ** Examples

  Kn(S=35,T=25,P=0, pHscale="T")



cleanEx()
nameEx("Ks")
### * Ks

flush(stderr()); flush(stdout())

### Name: Ks
### Title: Stability constant of hydrogen sulfate (mol/kg)
### Aliases: Ks
### Keywords: utilities

### ** Examples

  Ks(S=35,T=25,P=0, ks="d")



cleanEx()
nameEx("Ksi")
### * Ksi

flush(stderr()); flush(stdout())

### Name: Ksi
### Title: Dissociation constant of Si(OH)4
### Aliases: Ksi
### Keywords: utilities

### ** Examples

  Ksi(S=35, T=25, P=0, pHscale="T")



cleanEx()
nameEx("Kspa")
### * Kspa

flush(stderr()); flush(stdout())

### Name: Kspa
### Title: Solubility product of aragonite (mol/kg)
### Aliases: Kspa
### Keywords: utilities

### ** Examples

Kspa(S=35,T=25,P=0)



cleanEx()
nameEx("Kspc")
### * Kspc

flush(stderr()); flush(stdout())

### Name: Kspc
### Title: Solubility product of calcite (mol/kg)
### Aliases: Kspc
### Keywords: utilities

### ** Examples

	Kspc(S=35,T=25,P=0)



cleanEx()
nameEx("Kw")
### * Kw

flush(stderr()); flush(stdout())

### Name: Kw
### Title: Ion product of water (mol2/kg2)
### Aliases: Kw
### Keywords: utilities

### ** Examples

  Kw(S=35,T=25,P=0,pHscale="T")



cleanEx()
nameEx("Om")
### * Om

flush(stderr()); flush(stdout())

### Name: Om
### Title: Carbonate saturation state for magnesian calcites
### Aliases: Om
### Keywords: utilities

### ** Examples

Om(x=seq(0.01, 0.252, 0.01), flag=8, var1=8.2, var2=0.00234, 
  k1k2='x', kf='x', ks="d", pHscale="T", b="u74")



cleanEx()
nameEx("Pcorrect")
### * Pcorrect

flush(stderr()); flush(stdout())

### Name: Pcorrect
### Title: Pressure correction of equilibrium constants
### Aliases: Pcorrect
### Keywords: utilities

### ** Examples

k10 <- K1(T=25, P=0, S=35)
Pcorrect(Kvalue=k10, Ktype="K1", P=300, T=25, S=35, pHscale="T")



cleanEx()
nameEx("SIR")
### * SIR

flush(stderr()); flush(stdout())

### Name: sir
### Title: Parameters of the seawater carbonate system including the
###   substrate-inhibitor-ratio (SIR)
### Aliases: sir
### Keywords: utilities

### ** Examples

## With a couple of variables
sir(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## Using vectors as arguments
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
b <- c("l10", "l10", "l10")
sir(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P,
  Pt=Pt, Sit=Sit, kf=kf, k1k2=k1k2, pHscale=pHscale, b=b)

## Test with all flags 
flag <- c((1:15), (21:25))
var1 <- c(8.200000, 7.308171e-06, 7.308171e-06, 7.308171e-06, 7.308171e-06, 
	8.2, 8.2, 8.2, 8.2, 0.001646857, 0.001646857, 0.001646857, 0.0002822957, 
	0.0002822957, 0.00234, 258.2164, 258.2164, 258.2164, 258.2164, 258.2164 )
var2 <- c(7.308171e-06, 0.001646857, 0.0002822957, 0.00234, 0.001936461, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461, 0.0002822957, 
	0.00234, 0.001936461,  0.00234, 0.001936461, 0.001936461, 8.2, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461)
sir(flag=flag, var1=var1, var2=var2)

## Test using a data frame 
data(seacarb_test_P0)	#test data set for P=0 (surface)
tab <- seacarb_test_P0[14:19,]

## method 1 using the column numbers
sir(flag=tab[[1]], var1=tab[[2]], var2=tab[[3]], S=tab[[4]], T=tab[[5]], 
P=tab[[6]], Sit=tab[[8]], Pt=tab[[7]])

## method 2 using the column names
sir(flag=tab$flag, var1=tab$var1, var2=tab$var2, S=tab$S, T=tab$T, 
P=tab$P, Sit=tab$Sit, Pt=tab$Pt)



cleanEx()
nameEx("SIR_b")
### * SIR_b

flush(stderr()); flush(stdout())

### Name: sir_b
### Title: Parameters of the seawater carbonate system with boron addition,
###   including the substrate-inhibitor-ratio (SIR)
### Aliases: sir_b
### Keywords: utilities

### ** Examples

## With a couple of variables
sir_b(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, 
  Pt=0, Sit=0,pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", badd=0)




cleanEx()
nameEx("SIR_full")
### * SIR_full

flush(stderr()); flush(stdout())

### Name: sir_full
### Title: Parameters of the seawater carbonate system including the
###   substrate-inhibitor-ratio (SIR) - extension of carb
### Aliases: sir_full
### Keywords: utilities

### ** Examples

## With a couple of variables
sir_full(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="potential", NH4t=0, HSt=0)

## With a couple of variables and non-zero nutrient concentrations
sir_full(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=5e-6, Sit=2e-6,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="potential", NH4t=10e-6, HSt=0.1e-6)

## Using vectors as arguments
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
b <- c("l10", "l10", "l10")
gas <- c("potential", "potential", "potential")
NH4t <- c(0, 0, 0)
HSt <- c(0, 0, 0)
sir_full(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, Sit=Sit, 
  kf=kf, k1k2=k1k2, pHscale=pHscale, b=b, gas=gas, NH4t=NH4t, HSt=HSt)

## Test with all flags 
flag <- c((1:15), (21:25))
var1 <- c(8.200000, 7.308171e-06, 7.308171e-06, 7.308171e-06, 7.308171e-06, 
	8.2, 8.2, 8.2, 8.2, 0.001646857, 0.001646857, 0.001646857, 0.0002822957, 
	0.0002822957, 0.00234, 258.2164, 258.2164, 258.2164, 258.2164, 258.2164 )
var2 <- c(7.308171e-06, 0.001646857, 0.0002822957, 0.00234, 0.001936461, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461, 0.0002822957, 
	0.00234, 0.001936461,  0.00234, 0.001936461, 0.001936461, 8.2, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461)
sir_full(flag=flag, var1=var1, var2=var2)



cleanEx()
nameEx("amp")
### * amp

flush(stderr()); flush(stdout())

### Name: amp
### Title: pH value of the AMP buffer
### Aliases: amp
### Keywords: utilities

### ** Examples

	##Example from Dickson et al. (2007)
	amp(S=35,T=25)



cleanEx()
nameEx("bjerrum")
### * bjerrum

flush(stderr()); flush(stdout())

### Name: bjerrum
### Title: Bjerrum plot
### Aliases: bjerrum
### Keywords: utilities

### ** Examples

## Plot the bjerrum plot for the carbonate system using the default values
bjerrum(K1(),K2(),main="DIC speciation",lwd=2) 
abline(v=-log10(K1()),col="grey")
mtext(side=3,at=-log10(K1()),"pK1")
abline(v=-log10(K2()),col="grey")
mtext(side=3,at=-log10(K2()),"pK2")
legend("left",lty=1:3,lwd=2,legend=c(expression(CO[2]),expression(HCO[3]^"-"),
	expression(CO[3]^"2-")))

## Plot the bjerrum plot for phosphate using the default values
bjerrum(K1p(),K2p(),K3p(),main="phosphate speciation",lwd=2)
legend("left",lty=1:4,lwd=2,legend=c(expression(H[3]~PO[4]),
	expression(H[2]~PO[4]^"-"),
expression(HPO[4]^"2-"),expression(PO[4]^"3-")))

## Plot the bjerrum plot for the carbonate system using the values other 
##	than the default ones, showing the effect of temperature
bjerrum(K1(T=25,S=35),K2(T=25,S=35),conc=1.3,main="effect of temperature" )
bjerrum(K1(T=0,S=35),K2(T=0,S=35),conc=1.3,add=TRUE,col="red")
legend("left",lty=1,col=c("black","red"),legend=c("T=25 oC","T=0 oC"))
legend("right",lty=1:3,legend=c(expression(CO[2]),expression(HCO[3]^"-"),
	expression(CO[3]^"2-")))

## Plot the bjerrum plot for the carbonate system using the values other 
##	than the default ones, showing the effect of salinity
bjerrum(K1(T=25,S=35),K2(T=25,S=35),conc=1.3,main="effect of salinity" )
bjerrum(K1(T=25,S=5),K2(T=25,S=5),conc=1.3,add=TRUE,col="blue")
legend("left",lty=1,col=c("black","blue"),legend=c("S=35","S=5"))
legend("right",lty=1:3,legend=c(expression(CO[2]),expression(HCO[3]^"-"),
	expression(CO[3]^"2-")))

## Plot the bjerrum plot for the carbonate system using the values other 
##	than the default ones, showing the effect of pressure
bjerrum(K1(P=0),K2(P=0),conc=1.3,main="effect of pressure" )
bjerrum(K1(P=300),K2(P=300),conc=1.3,add=TRUE,col="green")
legend("left",lty=1,col=c("black","green"),legend=c("P=0","P=300"),title="atm")
legend("right",lty=1:3,legend=c(expression(CO[2]),expression(HCO[3]^"-"),
	expression(CO[3]^"2-")))



cleanEx()
nameEx("bor")
### * bor

flush(stderr()); flush(stdout())

### Name: bor
### Title: Total boron concentration (mol/kg)
### Aliases: bor
### Keywords: utilities

### ** Examples

bor(35, "l10")



cleanEx()
nameEx("buffer")
### * buffer

flush(stderr()); flush(stdout())

### Name: buffer
### Title: Buffer parameters of the seawater carbonate system
### Aliases: buffer
### Keywords: utilities

### ** Examples


## Computation with a couple of variables
buffer(flag=8, var1=8.2, var2=0.00234, S=35, T=25, Patm=1, P=0, Pt=0, 
	Sit=0, pHscale="T", kf="pf", k1k2="l", b="u74")

## Using vectors as arguments
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
b <- c("l10", "l10", "l10")
buffer(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, 
	Sit=Sit, kf=kf, k1k2=k1k2, pHscale=pHscale, b=b)


## Test for all flags 

flag <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 21, 22, 23, 24, 25)

var1 <- c(8.200000, 7.477544e-06, 7.477544e-06, 7.477544e-06, 7.477544e-06, 8.2, 
	8.2, 8.2, 8.2, 0.001685024, 0.001685024, 0.001685024,  0.0002888382, 
	0.0002888382, 0.002391252, 264.2008, 264.2008, 264.2008, 264.2008, 264.2008)

var2 <- c(7.477544e-06, 0.001685024, 0.0002888382, 0.002391252, 0.001981340, 
	0.001685024, 0.0002888382, 0.002391252, 0.001981340, 0.0002888382, 0.002391252,
	0.001981340,  0.002391252, 0.001981340, 0.001981340, 8.2, 0.001685024, 
	0.0002888382, 0.002391252, 0.001981340)

buffer(flag=flag, var1=var1, var2=var2)




cleanEx()
nameEx("buffergen")
### * buffergen

flush(stderr()); flush(stdout())

### Name: buffergen
### Title: Buffer factors of the seawater carbonate system as defined by
###   Hagens and Middelburg (2016)
### Aliases: buffergen
### Keywords: utilities

### ** Examples


## With a couple of variables
buffergen(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="potential", NH4t=0, HSt=0)

## With a couple of variables and non-zero nutrient concentrations
buffergen(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=5e-6, Sit=2e-6,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="potential", NH4t=10e-6, HSt=0.1e-6)

## Using vectors as arguments
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
b <- c("l10", "l10", "l10")
gas <- c("potential")
NH4t <- c(0, 0, 0)
HSt <- c(0, 0, 0)
buffergen(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, Sit=Sit, 
  kf=kf, k1k2=k1k2, pHscale=pHscale, b=b, gas=gas, NH4t=NH4t, HSt=HSt)

## Test with all flags 
flag <- c((1:15), (21:25))
var1 <- c(8.200000, 7.308171e-06, 7.308171e-06, 7.308171e-06, 7.308171e-06, 
	8.2, 8.2, 8.2, 8.2, 0.001646857, 0.001646857, 0.001646857, 0.0002822957, 
	0.0002822957, 0.00234, 258.2164, 258.2164, 258.2164, 258.2164, 258.2164 )
var2 <- c(7.308171e-06, 0.001646857, 0.0002822957, 0.00234, 0.001936461, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461, 0.0002822957, 
	0.00234, 0.001936461,  0.00234, 0.001936461, 0.001936461, 8.2, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461)
buffergen(flag=flag, var1=var1, var2=var2)



cleanEx()
nameEx("buffesm")
### * buffesm

flush(stderr()); flush(stdout())

### Name: buffesm
### Title: Buffer capacities of the seawater carbonate system from Egleston
###   et al. (2010), corrected and enhanced
### Aliases: buffesm
### Keywords: utilities

### ** Examples


## Computation with a couple of variables
buffesm(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Pt=0, 
	Sit=0, pHscale="T", kf="pf", k1k2="l", b="u74")

## Using vectors as arguments
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
b <- c("l10", "l10", "l10")
buffesm(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, 
	Sit=Sit, kf=kf, k1k2=k1k2, pHscale=pHscale, b=b)

## Test for all flags 
flag <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 21, 22, 23, 24, 25)

var1 <- c(8.200000, 7.477544e-06, 7.477544e-06, 7.477544e-06, 7.477544e-06, 8.2, 
	8.2, 8.2, 8.2, 0.001685024, 0.001685024, 0.001685024,  0.0002888382, 
	0.0002888382, 0.002391252, 264.2008, 264.2008, 264.2008, 264.2008, 264.2008)
var2 <- c(7.477544e-06, 0.001685024, 0.0002888382, 0.002391252, 0.001981340, 
	0.001685024, 0.0002888382, 0.002391252, 0.001981340, 0.0002888382, 0.002391252,
	0.001981340,  0.002391252, 0.001981340, 0.001981340, 8.2, 0.001685024, 
	0.0002888382, 0.002391252, 0.001981340)
buffesm(flag=flag, var1=var1, var2=var2)

## Compute 2 additional factors of interest (ratios of relative changes)
be <- buffesm(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, 
	Sit=Sit, kf=kf, k1k2=k1k2, pHscale=pHscale, b=b)
#     Ratio of gammaDIC/betaDIC = d ln [H+] / d ln pCO2
      Hfac <- (be$gammaDIC/be$betaDIC)                         #H+ factor
#     Ratio of gammaDIC/omegaDIC = d ln [CO32-] / d ln pCO2
      Satfac <- (be$gammaDIC/be$omegaDIC)                      #Saturation factor




cleanEx()
nameEx("carb")
### * carb

flush(stderr()); flush(stdout())

### Name: carb
### Title: Parameters of the seawater carbonate system
### Aliases: carb
### Keywords: utilities

### ** Examples


## With a couple of variables
carb(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## Using vectors as arguments
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
b <- c("l10", "l10", "l10")
carb(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P,
  Pt=Pt, Sit=Sit, kf=kf, k1k2=k1k2, pHscale=pHscale, b=b)

## Test with all flags 
flag <- c((1:15), (21:25))
var1 <- c(8.200000, 7.308171e-06, 7.308171e-06, 7.308171e-06, 7.308171e-06, 
	8.2, 8.2, 8.2, 8.2, 0.001646857, 0.001646857, 0.001646857, 0.0002822957, 
	0.0002822957, 0.00234, 258.2164, 258.2164, 258.2164, 258.2164, 258.2164 )
var2 <- c(7.308171e-06, 0.001646857, 0.0002822957, 0.00234, 0.001936461, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461, 0.0002822957, 
	0.00234, 0.001936461,  0.00234, 0.001936461, 0.001936461, 8.2, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461)
carb(flag=flag, var1=var1, var2=var2)

## Test using a data frame 
data(seacarb_test_P0)	#test data set for P=0 (surface)
tab <- seacarb_test_P0[14:19,]

## method 1 using the column numbers
carb(flag=tab[[1]], var1=tab[[2]], var2=tab[[3]], S=tab[[4]], T=tab[[5]], 
P=tab[[6]], Sit=tab[[8]], Pt=tab[[7]])

## method 2 using the column names
carb(flag=tab$flag, var1=tab$var1, var2=tab$var2, S=tab$S, T=tab$T, 
P=tab$P, Sit=tab$Sit, Pt=tab$Pt)




cleanEx()
nameEx("carbb")
### * carbb

flush(stderr()); flush(stdout())

### Name: carbb
### Title: Parameters of the seawater carbonate system with boron addition
### Aliases: carbb
### Keywords: utilities

### ** Examples


## With a couple of variables
carbb(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", badd=0)




cleanEx()
nameEx("carbfull")
### * carbfull

flush(stderr()); flush(stdout())

### Name: carbfull
### Title: Parameters of the seawater carbonate system - extension of carb
### Aliases: carbfull
### Keywords: utilities

### ** Examples


## With a couple of variables
carbfull(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="potential", NH4t=0, HSt=0)

## With a couple of variables and non-zero nutrient concentrations
carbfull(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=5e-6, Sit=2e-6,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="potential", NH4t=10e-6, HSt=0.1e-6)

## Using vectors as arguments
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
b <- c("l10", "l10", "l10")
gas <- c("potential", "potential", "potential")
NH4t <- c(0, 0, 0)
HSt <- c(0, 0, 0)
carbfull(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, Sit=Sit, 
  kf=kf, k1k2=k1k2, pHscale=pHscale, b=b, gas=gas, NH4t=NH4t, HSt=HSt)

## Test with all flags 
flag <- c((1:15), (21:25))
var1 <- c(8.200000, 7.308171e-06, 7.308171e-06, 7.308171e-06, 7.308171e-06, 
	8.2, 8.2, 8.2, 8.2, 0.001646857, 0.001646857, 0.001646857, 0.0002822957, 
	0.0002822957, 0.00234, 258.2164, 258.2164, 258.2164, 258.2164, 258.2164 )
var2 <- c(7.308171e-06, 0.001646857, 0.0002822957, 0.00234, 0.001936461, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461, 0.0002822957, 
	0.00234, 0.001936461,  0.00234, 0.001936461, 0.001936461, 8.2, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461)
carbfull(flag=flag, var1=var1, var2=var2)



cleanEx()
nameEx("d2p")
### * d2p

flush(stderr()); flush(stdout())

### Name: d2p
### Title: Converts depth in meters to pressure in dbar
### Aliases: d2p
### Keywords: utilities

### ** Examples

d2p(depth=7500, lat=30)  



cleanEx()
nameEx("derivnum")
### * derivnum

flush(stderr()); flush(stdout())

### Name: derivnum
### Title: Numerical derivatives of seawater carbonate system variables
### Aliases: derivnum
### Keywords: utilities

### ** Examples


## 1) For the input pair ALK and DIC (var1 and var2 when flag=15)
##    compute derivatives of all output varialbes 
##    with respect to DIC (i.e., var2)
derivnum(varid='var2', flag=15, var1=2300e-6, var2=2000e-6, 
         S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
         pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## 2) For the input pair pH and ALK (var1 and var2 when flag=8)
##    compute derivatives of all output variables 
##    with respect to [H+] concentration
derivnum(varid='var1', flag=8, var1=8.2, var2=0.00234, 
         S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
         pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## 3) Using vectors as arguments and compute derivatives of all output 
##    variables with respect to temperature
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
b <- c("u74", "u74", "u74")
derivnum(varid='T', flag=flag, var1=var1, var2=var2, S=S, T=T, P=P,
         Pt=Pt, Sit=Sit, kf=kf, k1k2=k1k2, pHscale=pHscale, b=b)

# For more examples of use of derivnum.R,
# consult the code of seacarb's errors routine.



cleanEx()
nameEx("eos2teos_chem")
### * eos2teos_chem

flush(stderr()); flush(stdout())

### Name: eos2teos_chem
### Title: Convert temperature and salinity from EOS-80 to TEOS-10
### Aliases: eos2teos_chem
### Keywords: utilities

### ** Examples

   # Calculate Conservative Temperature and Absolute Salinity of a sample with 
   # Practical salinity of 35 psu, in-situ temperature of 18 deg C,
   # at 0 dbar and total alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
   f <- eos2teos_chem(SP=35, T=18, P=0, TA=0.00234, DIC=0.00202)
   CT <- f$CT     # Conservative Temperature
   SA <- f$SA     # Absolute Salinity



cleanEx()
nameEx("eos2teos_geo")
### * eos2teos_geo

flush(stderr()); flush(stdout())

### Name: eos2teos_geo
### Title: Convert temperature and salinity from EOS-80 to TEOS-10
### Aliases: eos2teos_geo
### Keywords: utilities

### ** Examples

   # Calculate conservative temperature and absolute salinity of a sample with 
   # Practical salinity of 35 psu, in situ temperature of 18 deg C,
   # depth is 10 dbar and location is 188 degrees East and 4 degrees North.
   f <- eos2teos_geo(SP=35, T=18, P=10, long=188, lat=4)
   CT <- f$CT     # Conservative temperature
   SA <- f$SA     # Absolute salinity



cleanEx()
nameEx("errors")
### * errors

flush(stderr()); flush(stdout())

### Name: errors
### Title: Uncertainty propagation for computed marine carbonate system
###   variables
### Aliases: errors
### Keywords: utilities

### ** Examples


## 1) For the input pair ALK and DIC (var1 and var2 when flag=15),
## compute resulting uncertainty from given uncertainty on ALK and DIC (5 umol/kg)
## and default uncertainties in dissociation constants and total boron
## using the default method (Gaussian)
errors(flag=15, var1=2300e-6, var2=2000e-6, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0, 
       evar1=5e-6, evar2=5e-6, eS=0, eT=0, ePt=0, eSit=0, 
       pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")
## Typical output:
## H             pH          CO2           fCO2      pCO2      HCO3          ...
## 3.721614e-10  0.01796767  5.441869e-07  19.25338  19.31504  9.170116e-06  ...

## 2) Do the same as in one, but assign a 4% uncertainty to total boron
##    This uncertainty is the amount by which estimates from Lee et al (2010) and 
##    Uppstrom (1974) differ. The default for the latter is eBt=0.02.
errors(flag=15, var1=2300e-6, var2=2000e-6, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0, 
       evar1=5e-6, evar2=5e-6, eS=0, eT=0, ePt=0, eSit=0, eBt=0.04,
       pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## 3) For the input pair pH and ALK (var1 and var2 when flag=8)
## compute standard errors in output variables from errors in input variables, i.e., 
## for pH (0.005 pH units) and in ALK (5 umol/kg), along with
## errors in total dissolved inorganic phosphorus (0.1 umol/kg) and
## total dissolved inorganic silicon (2 umol/kg) concentrations, while
## assuming no uncertainty in dissociation constants & boron, using the Gaussian method:
errors(flag=8, var1=8.25, var2=2300e-6,  S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0, 
       evar1=0.005, evar2=5e-6, eS=0, eT=0, ePt=0.1, eSit=2, epK=0, eBt=0,
       method="ga", pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## 4) For the input pair pCO2 and pH (var1 and var2 when flag=21)
## compute standard errors in output variables from errors in input variables, i.e., 
## for pCO2 (2 uatm) and pH (0.005 pH units), with no uncertainties in Pt and Sit
## nor in the dissociation constants BUT a perfect anticorrelation between pCO2 and pH,
## (the input pair) using the Method of moments:
errors(flag=21, var1=400, var2=8.1,  S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0, 
       evar1=2, evar2=0.005, eS=0, eT=0, ePt=0.0, eSit=0, epK=0, eBt=0, 
       method="mo", r=-1.0, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## 5) Use vectors as arguments and compute errors on all output variables
## using Monte Carlo method taking into account input errors on pH, ALK, DIC
## and dissociation constants (pKx)
flag <- c(8, 15)
var1 <- c(8.2, 0.002394, 8.25)
var2 <- c(0.002343955, 0.002017)
S <- c(35, 35)
T <- c(25, 25)
P <- 0
Pt <- 0
Sit <- 0
evar1 <- c(0.005, 2e-6)
evar2 <- c(2e-6, 2e-6)
epKx <- c(0.002, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01)
eBtx = 0.01
method <- "mc"
kf <- "pf"
k1k2 <- "l"
pHscale <- "T"
b <- "u74"
## NOTE that the following is executable but enclosed in "donttest" 
## because it takes too long to run when submiting to CRAN
## and is therefore rejected



cleanEx()
nameEx("f2pCO2")
### * f2pCO2

flush(stderr()); flush(stdout())

### Name: f2pCO2
### Title: Converts the CO2 fugacity to CO2 partial pressure
### Aliases: f2pCO2
### Keywords: utilities

### ** Examples

f2pCO2(T=25, Patm=1.0, P=0, fCO2=380)  



cleanEx()
nameEx("fCO2insi")
### * fCO2insi

flush(stderr()); flush(stdout())

### Name: fCO2insi
### Title: fCO2 at in situ temperature
### Aliases: fCO2insi
### Keywords: utilities

### ** Examples

fCO2insi(fCO2lab = 400, SST = 15, Tlab = 16)



cleanEx()
nameEx("kconv")
### * kconv

flush(stderr()); flush(stdout())

### Name: kconv
### Title: Conversion factors to change the pH scale of dissociation
###   constants
### Aliases: kconv
### Keywords: utilities

### ** Examples

##To convert dissociation constants from the total scale to the free scale
## (at salinity=35, temperature=25oC and atmospheric pressure):
kconv(35,25,0)
conv <- kconv()
 c(K1_total=K1(),K1_SWS=K1()*conv$ktotal2SWS,K1_free=K1()*conv$ktotal2free)  



cleanEx()
nameEx("oa")
### * oa

flush(stderr()); flush(stdout())

### Name: oa
### Title: Perturbation of the seawater carbonate system
### Aliases: oa
### Keywords: utilities

### ** Examples

oa(flag=24, var1=384, var2=2325e-6, pCO2s=1e6, pCO2f=793, S=34.3, T=16, 
	P=0, pHscale="T", kf="pf", k1k2="l", ks="d", plot=TRUE, b="u74")



cleanEx()
nameEx("p2d")
### * p2d

flush(stderr()); flush(stdout())

### Name: p2d
### Title: Converts pressure in dbar to depth in meters
### Aliases: p2d
### Keywords: utilities

### ** Examples

p2d(pressure=7686, lat=30)  



cleanEx()
nameEx("p2fCO2")
### * p2fCO2

flush(stderr()); flush(stdout())

### Name: p2fCO2
### Title: Converts pCO2 (partial pressure in CO2) into fCO2 (fugacity of
###   CO2)
### Aliases: p2fCO2
### Keywords: utilities

### ** Examples

p2fCO2(T=25, Patm=0.97, P=0, pCO2=380)  



cleanEx()
nameEx("p2xCO2")
### * p2xCO2

flush(stderr()); flush(stdout())

### Name: p2xCO2
### Title: Converts partial pressure of CO2 to mole fraction of CO2
### Aliases: p2xCO2
### Keywords: utilities

### ** Examples

## Convert atmospheric pressure from mbar to atm
   Patm_mbar = 1052                 # in millibar
   Patm      = Patm_mbar / 1013.25  # in atm
## Compute xCO2 from pCO2
   pCO2 = 380
   xCO2 = p2xCO2(T=25, S=35, Patm=Patm, pCO2=pCO2)
   print(xCO2)
## The result is 377.1546 ppm 



cleanEx()
nameEx("pCa")
### * pCa

flush(stderr()); flush(stdout())

### Name: pCa
### Title: pCa
### Aliases: pCa
### Keywords: utilities

### ** Examples

pCa(flag=15, var1=2302e-6, var2=2050e-6, Ca=0.01028, S=35, T=20, P=0, 
	Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74") # with normal Ca concentration
pCa(flag=15, var1=2302e-6, var2=2050e-6, Ca=0.01028/2, S=35, T=20, P=0, 
	Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74") # with 0.5 * Ca concentration



cleanEx()
nameEx("pH")
### * pH

flush(stderr()); flush(stdout())

### Name: pH
### Title: Potentiometric pH
### Aliases: pH
### Keywords: utilities

### ** Examples

	##Example from Dickson et al. (2007)
	pH(Ex=-67,Etris=-72.4,S=35,T=25)



cleanEx()
nameEx("pHconv")
### * pHconv

flush(stderr()); flush(stdout())

### Name: pHconv
### Title: Conversion of pH
### Aliases: pHconv
### Keywords: utilities

### ** Examples

##To convert pH=8.10 from the seawater scale to the total scale
##at salinity=35, temperature=25oC and atmospheric pressure:

  pHc <- pHconv(flag=1, pH=8.10, S=35, T=25, P=0, ks="d")

##note that pHc is the value of the pH converted in total scale

## By using vectors
## to convert the pH values : 8, 8.05, 8.10, 8.15, 8.20 
## from the free to the total scale

pH <- c(8, 8.05, 8.10, 8.15, 8.20)
pHc <- pHconv(flag=2, pH=pH, S=35, T=25, P=0, ks="d")

## note that pHc is a vector containing the value of the pH converted 
## to the total scale



cleanEx()
nameEx("pHinsi")
### * pHinsi

flush(stderr()); flush(stdout())

### Name: pHinsi
### Title: pH at in situ temperature and pressure
### Aliases: pHinsi
### Keywords: utilities

### ** Examples

  pHinsi(pH=8.2,ALK=2.4e-3,Tinsi=25,Tlab=25,Pinsi=200,S=35,Pt=0,Sit=0)



cleanEx()
nameEx("pHslope")
### * pHslope

flush(stderr()); flush(stdout())

### Name: pHslope
### Title: Slope of the calibration curve of a pH electrode
### Aliases: pHslope
### Keywords: utilities

### ** Examples

	##Example from Dickson et al. (2007)
	pHslope(Etris=-72.4,Eamp=4.9,S=35,T=25)



cleanEx()
nameEx("pHspec")
### * pHspec

flush(stderr()); flush(stdout())

### Name: pHspec
### Title: Calculates pHT from results of spectrophotometric measurements
### Aliases: pHspec
### Keywords: utilities

### ** Examples

	##Example should give test value pHT = 7.6713
	 pHspec(S=35, T=25, R=1, d="mCP", k="m18", warn="y")



cleanEx()
nameEx("pTA")
### * pTA

flush(stderr()); flush(stdout())

### Name: pTA
### Title: pTA
### Aliases: pTA
### Keywords: utilities

### ** Examples

pTA(flag=24, sys=0, var1=384, var2=2302e-6, pCO2a=384, co3=260e-6, 
	hco3=1000e-6, S=34.3, T=16, P=0, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

pTA(flag=24, sys=1, var1=384, var2=2302e-6, pCO2a=384, co3=260e-6, 
	hco3=1000e-6, S=34.3, T=16, P=0, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")



cleanEx()
nameEx("pgas")
### * pgas

flush(stderr()); flush(stdout())

### Name: pgas
### Title: pgas
### Aliases: pgas
### Keywords: utilities

### ** Examples

pgas(flag=15, var1=2302e-6, var2=2050e-6, pCO2g=750, S=35, T=20, P=0, 
	Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74") 



cleanEx()
nameEx("pmix")
### * pmix

flush(stderr()); flush(stdout())

### Name: pmix
### Title: pmix
### Aliases: pmix
### Keywords: utilities

### ** Examples

pmix(flag=24, var1=384, var2=2302e-6, pCO2s=1e6, wf=0.003, S=34.3, 
	T=16, P=0, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")



cleanEx()
nameEx("ppH")
### * ppH

flush(stderr()); flush(stdout())

### Name: ppH
### Title: ppH
### Aliases: ppH
### Keywords: utilities

### ** Examples

ppH(flag=24, sys=0, var1=384, var2=2302e-6, pCO2a=384, vol=-12e-3, 
	N=0.01, S=34.3, T=16, P=0, pHscale="T", kf="pf", k1k2="l", ks="d")

ppH(flag=24, sys=1, var1=384, var2=2302e-6, pCO2a=384, vol=-12e-3, 
	N=0.01, S=34.3, T=16, P=0, pHscale="T", kf="pf", k1k2="l", ks="d")



cleanEx()
nameEx("psi")
### * psi

flush(stderr()); flush(stdout())

### Name: psi
### Title: Molar ratio of CO2 released vs CaCO3 precipitated
### Aliases: psi
### Keywords: utilities

### ** Examples

## Calculation using the numerical example given in Frankignoulle et al. (1994)
psi(flag=24, var1=350, var2=2400e-6, S=35, T=25, P=0, Pt=0, 
	Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d")



cleanEx()
nameEx("rho")
### * rho

flush(stderr()); flush(stdout())

### Name: rho
### Title: Density of seawater (kg/m3)
### Aliases: rho
### Keywords: utilities

### ** Examples

	rho(35,25,0)



cleanEx()
nameEx("sa2sp_chem")
### * sa2sp_chem

flush(stderr()); flush(stdout())

### Name: sa2sp_chem
### Title: From absolute to practical salinity
### Aliases: sa2sp_chem
### Keywords: utilities

### ** Examples

   # Calculate the practical salinity of a sample with Absolute Salinity of 35 g/kg,
   # Total Alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
   SP <- sa2sp_chem(SA=35, TA=0.00234, DIC=0.00202)



cleanEx()
nameEx("sa2sp_geo")
### * sa2sp_geo

flush(stderr()); flush(stdout())

### Name: sa2sp_geo
### Title: From absolute to practical salinity
### Aliases: sa2sp_geo
### Keywords: utilities

### ** Examples

   # Calculate the practical salinity of a sample whose absolute Salinity is 35,
   # depth is 10 dbar and location is 188 degrees East and 4 degrees North.
   SP <- sa2sp_geo(35, 10, 188, 4)



cleanEx()
nameEx("sf_calc")
### * sf_calc

flush(stderr()); flush(stdout())

### Name: sf_calc
### Title: Calculation of calibrated pH for seaFET sensor
### Aliases: sf_calc
### Keywords: utilities

### ** Examples

sf_calc(calEint=0.0865, calEext= -0.93, E0int25 =-0.39, 
  E0ext25=-1.46, calT=16.2, calSal=35.6)

## Using the test file seaFET
sf_calc(calEint=seaFET$Eint, calEext=seaFET$Eext, 
  E0int25=seaFET$E0int25, E0ext25=seaFET$E0ext25, 
  calT=seaFET$Temperature, calSal=seaFET$Salinity)



cleanEx()
nameEx("sf_calib")
### * sf_calib

flush(stderr()); flush(stdout())

### Name: sf_calib
### Title: Calibration coefficients for seaFET sensor
### Aliases: sf_calib
### Keywords: utilities

### ** Examples

sf_calib(calEint=0.0865, calEext=-0.93, calpH=8.132, calT=16.2, calSal=35.6)

## Using the test file seaFET
sf_calib(calEint=seaFET$Eint, calEext=seaFET$Eext, 
  calpH=seaFET$pHspectro, calT=seaFET$Temperature, 
  calSal=seaFET$Salinity)



cleanEx()
nameEx("sp2sa_chem")
### * sp2sa_chem

flush(stderr()); flush(stdout())

### Name: sp2sa_chem
### Title: From Practical to absolute salinity
### Aliases: sp2sa_chem
### Keywords: utilities

### ** Examples

   # Calculate the absolute salinity of a sample with practical Salinity of 35,
   # Total Alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
   SA <- sp2sa_chem(SP=35, TA=0.00234, DIC=0.00202)



cleanEx()
nameEx("sp2sa_geo")
### * sp2sa_geo

flush(stderr()); flush(stdout())

### Name: sp2sa_geo
### Title: From practical to absolute salinity
### Aliases: sp2sa_geo
### Keywords: utilities

### ** Examples

   # Calculate the absolute salinity of a sample whose practical Salinity is 35,
   # depth is 10 dbar and location is 188 degrees East and 4 degrees North.
   SA <- sp2sa_geo(35, 10, 188, 4)     # 34.711778344814114



cleanEx()
nameEx("speciation")
### * speciation

flush(stderr()); flush(stdout())

### Name: speciation
### Title: ionic forms as a function of pH
### Aliases: speciation
### Keywords: utilities

### ** Examples

## Speciation of divalent species; example to estimate the various ionic forms
## of dissolved inorganic carbon  (DIC = 0.0021 mol/kg) at a salinity of 35,
## a temperature of 25oC and an hydrostatic pressure of 0:
spec <- speciation (K1(35, 25, 0), K2(35, 25, 0), pH=8, conc=0.0021)
## where (spec\$C1=[CO2], spec\$C2=[HCO3-], spec\$C3=[CO3--])

## Speciation of trivalent species (e.g.,  H3PO4, H2PO4-, HPO4--, PO4---)
speciation(K1p(), K2p(), K3p(), conc=0.001)

## Effect of temperature on pCO2 - Figure 1.4.18 of Zeebe and Wolf-Gladrow (2001)
Tseq <- seq(0, 30, by=0.5)
pHseq <- carb(flag=15, var1=2300e-6, var2=1900e-6, S=35, T=Tseq, P=0)$pH
CO2  <- speciation(K1(T=Tseq), K2(T=Tseq), conc=1900, pH=pHseq)$C1
pCO2 <- CO2/K0(T=Tseq)
plot(Tseq, pCO2, xlab="Temperature (oC)", ylab="pCO2 (uatm)", type="l", 
	main="effect of temperature on pCO2")
legend("topleft", c(expression(sum(CO[2])==1900~umol~kg^"-1"), 
	expression(TA==2300~umol~kg^"-1")))



cleanEx()
nameEx("teos2eos_chem")
### * teos2eos_chem

flush(stderr()); flush(stdout())

### Name: teos2eos_chem
### Title: Convert temperature and salinity from TEOS-10 to EOS-80
### Aliases: teos2eos_chem
### Keywords: utilities

### ** Examples

   # Calculate in situ temperature and practical salinity of a sample with 
   # Absolute salinity of 35 g/kg, Conservative temperature of 18 deg C,
   # at 0 dbar and Total alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
   f <- teos2eos_chem(SA=35, CT=18, P=0, TA=0.00234, DIC=0.00202)
   T <- f$T     # insitu temperature
   SP <- f$SP     # Practical salinity



cleanEx()
nameEx("teos2eos_geo")
### * teos2eos_geo

flush(stderr()); flush(stdout())

### Name: teos2eos_geo
### Title: Convert temperature and salinity from TEOS-10 to EOS-80
### Aliases: teos2eos_geo
### Keywords: utilities

### ** Examples

   # Calculate in situ temperature and practical salinity of a sample with 
   # Absolute salinity of 35 g/kg, conservative temperature of 18 deg C,
   # depth is 10 dbar and location is 188 degrees East and 4 degrees North.
   f <- teos2eos_geo(SA=35, CT=18, P=10, long=188, lat=4)
   T <- f$T     # in situ temperature
   SP <- f$SP     # Practical salinity



cleanEx()
nameEx("theta")
### * theta

flush(stderr()); flush(stdout())

### Name: theta
### Title: Potential temperature of seawater
### Aliases: theta
### Keywords: utilities

### ** Examples

   #Calculate the potential temperature for a sample at 1000 db referenced to the surface
   theta <- theta(S=35, T=25, P=100, Pref=0)



cleanEx()
nameEx("tris")
### * tris

flush(stderr()); flush(stdout())

### Name: tris
### Title: pH of TRIS buffer
### Aliases: tris
### Keywords: utilities

### ** Examples

	##Example from Mueller et al. (2018), should give test value pHT = 8.0703
	 tris(S=20,T=25,b=0.04,k="m18")



cleanEx()
nameEx("vapress")
### * vapress

flush(stderr()); flush(stdout())

### Name: vapress
### Title: Computes vapor pressure of seawater
### Aliases: vapress
### Keywords: utilities

### ** Examples

   pH20 <- vapress(S=35, T=25, form="d2007")



cleanEx()
nameEx("x2pCO2")
### * x2pCO2

flush(stderr()); flush(stdout())

### Name: x2pCO2
### Title: Converts mole fraction to partial pressure of CO2
### Aliases: x2pCO2
### Keywords: utilities

### ** Examples

## Atmospheric pressure is rarely equal to 1 atm exactly 
## Over the Southern Ocean Patm=0.97 is more realistic
   pCO2_socn <- x2pCO2(S=35, T=0, Patm=0.97, xCO2=400.0)
   print(pCO2_socn)
## The result (385.6322 uatm) is 12 uatm less than if it was wrongly assumed that Patm=1.0

## Show effect of temperature on pCO2 computed from xCO2, and on resulting variables from "carb"
S <- 35
ALK <- 2300e-6
T <- seq(0,30,5)
xCO2 <- 400
pCO2 <- x2pCO2(S=35, T=T, Patm=1, xCO2=400)
results <- carb(flag=24, var1=pCO2, var2=ALK, S=S, T=T, P=0, Pt=0, Sit=0, 
  pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")
print(results)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
