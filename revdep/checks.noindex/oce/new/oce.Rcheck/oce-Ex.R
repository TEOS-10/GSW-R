pkgname <- "oce"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('oce')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CTD_BCD2014666_008_1_DN.ODF.gz")
### * CTD_BCD2014666_008_1_DN.ODF.gz

flush(stderr()); flush(stdout())

### Name: CTD_BCD2014666_008_1_DN.ODF.gz
### Title: Sample ctd dataset in odf format
### Aliases: CTD_BCD2014666_008_1_DN.ODF.gz

### ** Examples

ctd <- read.ctd(system.file("extdata", "CTD_BCD2014666_008_1_DN.ODF.gz", package="oce"))
summary(ctd)
plot(ctd)




cleanEx()
nameEx("GMTOffsetFromTz")
### * GMTOffsetFromTz

flush(stderr()); flush(stdout())

### Name: GMTOffsetFromTz
### Title: Determine time offset from timezone
### Aliases: GMTOffsetFromTz

### ** Examples

library(oce)
cat("Atlantic Standard Time is ", GMTOffsetFromTz("AST"), "hours after UTC")



cleanEx()
nameEx("T68fromT90")
### * T68fromT90

flush(stderr()); flush(stdout())

### Name: T68fromT90
### Title: Convert from ITS-90 to IPTS-68 temperature
### Aliases: T68fromT90

### ** Examples

library(oce)
T68 <- seq(3, 20, 1)
T90 <- T90fromT68(T68)
sqrt(mean((T68-T90)^2))




cleanEx()
nameEx("T90fromT48")
### * T90fromT48

flush(stderr()); flush(stdout())

### Name: T90fromT48
### Title: Convert from ITS-48 to ITS-90 temperature
### Aliases: T90fromT48

### ** Examples

library(oce)
T68 <- seq(3, 20, 1)
T90 <- T90fromT68(T68)
sqrt(mean((T68-T90)^2))




cleanEx()
nameEx("T90fromT68")
### * T90fromT68

flush(stderr()); flush(stdout())

### Name: T90fromT68
### Title: Convert from IPTS-68 to ITS-90 temperature
### Aliases: T90fromT68

### ** Examples

library(oce)
T68 <- seq(3, 20, 1)
T90 <- T90fromT68(T68)
sqrt(mean((T68-T90)^2))




cleanEx()
nameEx("ad2cpHeaderValue")
### * ad2cpHeaderValue

flush(stderr()); flush(stdout())

### Name: ad2cpHeaderValue
### Title: Decode an item from a Nortek AD2CP file header
### Aliases: ad2cpHeaderValue

### ** Examples

## Not run: 
##D d <- read.oce("a.ad2cp")
##D # The examples start with the line in x[["text"]][[1]]; note that in the second
##D # example, it would be insuficient to use a key of "BEAMCFGLIST", because that will
##D # yield 4 lines, and the function is not designed to handle that.
##D 
##D # ID,STR=\"Signature1000\",SN=123456
##D type <- ad2cpHeaderValue(d, "ID", "STR", numeric=FALSE)
##D serialNumber <- ad2cpHeaderValue(d, "ID", "SN")
##D 
##D # BEAMCFGLIST,BEAM=1,THETA=25.00,PHI=0.00,FREQ=1000,BW=25,BRD=1,HWBEAM=1,ZNOM=60.00
##D beam1Angle <- ad2cpHeaderValue(d, "BEAMCFGLIST,BEAM=1", "THETA")
##D frequency <- ad2cpHeaderValue(d, "BEAMCFGLIST,BEAM=1", "FREQ", default=NA)
## End(Not run)




cleanEx()
nameEx("addSpine")
### * addSpine

flush(stderr()); flush(stdout())

### Name: addSpine
### Title: Add a spine to a section object
### Aliases: addSpine

### ** Examples

library(oce)
data(section)
sectionWest <- subset(section, longitude < -60)
spine <- list(longitude=c(-74.5, -69.2, -55), latitude=c(38.6, 36.25, 36.25))
sectionWithSpine <- addSpine(sectionWest, spine)
plot(sectionWithSpine, which="map")
plot(sectionWithSpine, xtype="distance", which="temperature")
plot(sectionWithSpine, xtype="spine", which="temperature")




cleanEx()
nameEx("adp")
### * adp

flush(stderr()); flush(stdout())

### Name: adp
### Title: Sample adp (acoustic-doppler profiler) dataset
### Aliases: adp

### ** Examples






cleanEx()
nameEx("adpEnsembleAverage")
### * adpEnsembleAverage

flush(stderr()); flush(stdout())

### Name: adpEnsembleAverage
### Title: Ensemble Average an ADP Object in Time
### Aliases: adpEnsembleAverage

### ** Examples

library(oce)
data(adp)
adpAvg <- adpEnsembleAverage(adp, n=2)
plot(adpAvg)




cleanEx()
nameEx("adp_rdi.000")
### * adp_rdi.000

flush(stderr()); flush(stdout())

### Name: adp_rdi.000
### Title: Sample adp (acoustic-doppler profiler) file in RDI format
### Aliases: adp_rdi.000

### ** Examples

## Not run: 
##D read.oce(system.file("extdata", "adp_rdi.000", package="oce"))
## End(Not run)




cleanEx()
nameEx("adv-class")
### * adv-class

flush(stderr()); flush(stdout())

### Name: adv-class
### Title: Class to Store adv Data
### Aliases: adv-class

### ** Examples

data(adv)
adv[["v"]] <- 0.001 + adv[["v"]] # add 1mm/s to all velocity components




cleanEx()
nameEx("adv")
### * adv

flush(stderr()); flush(stdout())

### Name: adv
### Title: Sample adv (acoustic-doppler velocimeter) dataset
### Aliases: adv

### ** Examples





cleanEx()
nameEx("airRho")
### * airRho

flush(stderr()); flush(stdout())

### Name: airRho
### Title: Air density
### Aliases: airRho

### ** Examples

degC <- seq(0,30,length.out=100)
p <- seq(98,102,length.out=100) * 1e3
contour(x=degC, y=p, z=outer(degC,p,airRho), labcex=1)



cleanEx()
nameEx("amsr")
### * amsr

flush(stderr()); flush(stdout())

### Name: amsr
### Title: An amsr dataset for waters near Nova Scotia
### Aliases: amsr

### ** Examples

library(oce)
data(coastlineWorld)
data(amsr)
plot(amsr, "SST")
lines(coastlineWorld[["longitude"]], coastlineWorld[["latitude"]])




cleanEx()
nameEx("angle2hms")
### * angle2hms

flush(stderr()); flush(stdout())

### Name: angle2hms
### Title: Convert astronomical angle in degrees to hours, minutes and
###   seconds
### Aliases: angle2hms

### ** Examples

# A randomly-chosen example on page 99 of Meeus (1991).
angle2hms(177.74208) # string component 11h50m58s.10




cleanEx()
nameEx("angleRemap")
### * angleRemap

flush(stderr()); flush(stdout())

### Name: angleRemap
### Title: Convert angles from 0:360 to -180:180
### Aliases: angleRemap

### ** Examples


library(oce)
## fake some heading data that lie near due-north (0 degrees)
n <- 20
heading <- 360 + rnorm(n, sd=10)
heading <- ifelse(heading > 360, heading - 360, heading)
x <- 1:n
plot(x, heading, ylim=c(-10, 360), type='l', col='lightgray', lwd=10)
lines(x, angleRemap(heading))



cleanEx()
nameEx("approx3d")
### * approx3d

flush(stderr()); flush(stdout())

### Name: approx3d
### Title: Trilinear interpolation in a 3D array
### Aliases: approx3d

### ** Examples

## set up a grid
library(oce)
n <- 5
x <- seq(0, 1, length.out=n)
y <- seq(0, 1, length.out=n)
z <- seq(0, 1, length.out=n)
f <- array(1:n^3, dim=c(length(x), length(y), length(z)))
## interpolate along a diagonal line
m <- 100
xout <- seq(0, 1, length.out=m)
yout <- seq(0, 1, length.out=m)
zout <- seq(0, 1, length.out=m)
approx <- approx3d(x, y, z, f, xout, yout, zout)
## graph the results
plot(xout, approx, type='l')
points(xout[1], f[1, 1, 1])
points(xout[m], f[n,n,n])



cleanEx()
nameEx("argo")
### * argo

flush(stderr()); flush(stdout())

### Name: argo
### Title: ARGO float dataset
### Aliases: argo

### ** Examples

library(oce)
data(argo)
summary(argo)
data(coastlineWorld)
plot(argo, which="trajectory")




cleanEx()
nameEx("argoGrid")
### * argoGrid

flush(stderr()); flush(stdout())

### Name: argoGrid
### Title: Grid Argo float data
### Aliases: argoGrid

### ** Examples

library(oce)
data(argo)
g <- argoGrid(argo, p=seq(0, 100, 1))
par(mfrow=c(2,1))
t <- g[["time"]]
z <- -g[["pressure"]][,1]
## Set zlim because of spurious temperatures.
imagep(t, z, t(g[['temperature']]), ylim=c(-100,0), zlim=c(0,20))
imagep(t, z, t(g[['salinity']]), ylim=c(-100,0))




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("argoJuldToTime")
### * argoJuldToTime

flush(stderr()); flush(stdout())

### Name: argoJuldToTime
### Title: Convert Argo Julian Day (juld) to time
### Aliases: argoJuldToTime

### ** Examples

argoJuldToTime(25749)




cleanEx()
nameEx("as.adp")
### * as.adp

flush(stderr()); flush(stdout())

### Name: as.adp
### Title: Create an ADP Object
### Aliases: as.adp

### ** Examples

data(adp)
t <- adp[["time"]]
d <- adp[["distance"]]
v <- adp[["v"]]
a <- as.adp(time=t, distance=d, v=v)




cleanEx()
nameEx("as.cm")
### * as.cm

flush(stderr()); flush(stdout())

### Name: as.cm
### Title: Coerce data into a CM object
### Aliases: as.cm

### ** Examples

library(oce)
# Example 1: creation from scratch
t <- Sys.time() + 0:50
u <- sin(2*pi*0:50/5) + rnorm(51)
v <- cos(2*pi*0:50/5) + rnorm(51)
p <- 100 + rnorm(51)
summary(as.cm(t, u, v, p))

# Example 2: creation from an adv object
data(adv)
summary(as.cm(adv))




cleanEx()
nameEx("as.ctd")
### * as.ctd

flush(stderr()); flush(stdout())

### Name: as.ctd
### Title: Coerce data into CTD object
### Aliases: as.ctd

### ** Examples

library(oce)
## 1. fake data, with default units
pressure <- 1:50
temperature <- 10 - tanh((pressure - 20) / 5) + 0.02*rnorm(50)
salinity <- 34 + 0.5*tanh((pressure - 20) / 5) + 0.01*rnorm(50)
ctd <- as.ctd(salinity, temperature, pressure)
# Add a new column
fluo <- 5 * exp(-pressure / 20)
ctd <- oceSetData(ctd, name="fluorescence", value=fluo,
                  unit=list(unit=expression(mg/m^3), scale=""))
summary(ctd)

## 2. fake data, with supplied units (which are the defaults, actually)
ctd <- as.ctd(salinity, temperature, pressure,
    units=list(salinity=list(unit=expression(), scale="PSS-78"),
    temperature=list(unit=expression(degree*C), scale="ITS-90"),
    pressure=list(unit=expression(dbar), scale="")))




cleanEx()
nameEx("as.sealevel")
### * as.sealevel

flush(stderr()); flush(stdout())

### Name: as.sealevel
### Title: Coerce Data Into a Sealevel Object
### Aliases: as.sealevel

### ** Examples

library(oce)

# Construct a year of M2 tide, starting at the default time
# 0000-01-01T00:00:00.
h <- seq(0, 24*365)
elevation <- 2.0 * sin(2*pi*h/12.4172)
sl <- as.sealevel(elevation)
summary(sl)

# As above, but start at the Y2K time.
time <- as.POSIXct("2000-01-01") + h * 3600
sl <- as.sealevel(elevation, time)
summary(sl)



cleanEx()
nameEx("as.section")
### * as.section

flush(stderr()); flush(stdout())

### Name: as.section
### Title: Create a Section
### Aliases: as.section

### ** Examples

library(oce)
data(ctd)
## vector of names of CTD objects
fake <- ctd
fake[["temperature"]] <- ctd[["temperature"]] + 0.5
fake[["salinity"]] <- ctd[["salinity"]] + 0.1
fake[["longitude"]] <- ctd[["longitude"]] + 0.01
fake[["station"]] <- "fake"
sec1 <- as.section(c("ctd", "fake"))
summary(sec1)
## vector of CTD objects
ctds <- vector("list", 2)
ctds[[1]] <- ctd
ctds[[2]] <- fake
sec2 <- as.section(ctds)
summary(sec2)
## argo data (a subset)
data(argo)
sec3 <- as.section(subset(argo, profile<5))
summary(sec3)




cleanEx()
nameEx("as.tidem")
### * as.tidem

flush(stderr()); flush(stdout())

### Name: as.tidem
### Title: Create tidem object from fitted harmonic data
### Aliases: as.tidem

### ** Examples

# Simulate a tide table with output from tidem().
data(sealevelTuktoyaktuk)
# 'm0' is model fitted by tidem()
m0 <- tidem(sealevelTuktoyaktuk)
p0 <- predict(m0, sealevelTuktoyaktuk[["time"]])
m1 <- as.tidem(mean(sealevelTuktoyaktuk[["time"]]), sealevelTuktoyaktuk[["latitude"]],
               m0[["name"]], m0[["amplitude"]], m0[["phase"]])
# Test agreement with tidem() result, by comparing predicted sealevels.
p1 <- predict(m1, sealevelTuktoyaktuk[["time"]])
expect_lt(max(abs(p1 - p0), na.rm=TRUE), 1e-10)
# Simplified harmonic model, using large constituents
# > m0[["name"]][which(m[["amplitude"]]>0.05)]
# [1] "Z0"  "MM"  "MSF" "O1"  "K1"  "OO1" "N2"  "M2"  "S2"
h <- "
name  amplitude      phase
  Z0 1.98061875   0.000000
  MM 0.21213065 263.344739
 MSF 0.15605629 133.795004
  O1 0.07641438  74.233130
  K1 0.13473817  81.093134
 OO1 0.05309911 235.749693
  N2 0.08377108  44.521462
  M2 0.49041340  77.703594
  S2 0.22023705 137.475767"
coef <- read.table(text=h, header=TRUE)
m2 <- as.tidem(mean(sealevelTuktoyaktuk[["time"]]),
               sealevelTuktoyaktuk[["latitude"]],
               coef$name, coef$amplitude, coef$phase)
p2 <- predict(m2, sealevelTuktoyaktuk[["time"]])
expect_lt(max(abs(p2 - p0), na.rm=TRUE), 1)
par(mfrow=c(3, 1))
oce.plot.ts(sealevelTuktoyaktuk[["time"]], p0)
ylim <- par("usr")[3:4] # to match scales in other panels
oce.plot.ts(sealevelTuktoyaktuk[["time"]], p1, ylim=ylim)
oce.plot.ts(sealevelTuktoyaktuk[["time"]], p2, ylim=ylim)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("as.unit")
### * as.unit

flush(stderr()); flush(stdout())

### Name: as.unit
### Title: Convert a String to a Unit
### Aliases: as.unit

### ** Examples

as.unit("DBAR")
as.unit("IPTS-68")
as.unit("ITS-90")
as.unit("PSS-78")
as.unit("UMOL/KG")



cleanEx()
nameEx("as.windrose")
### * as.windrose

flush(stderr()); flush(stdout())

### Name: as.windrose
### Title: Create a Windrose Object
### Aliases: as.windrose

### ** Examples

library(oce)
xcomp <- rnorm(360) + 1
ycomp <- rnorm(360)
wr <- as.windrose(xcomp, ycomp)
summary(wr)
plot(wr)




cleanEx()
nameEx("bcdToInteger")
### * bcdToInteger

flush(stderr()); flush(stdout())

### Name: bcdToInteger
### Title: Decode BCD to integer
### Aliases: bcdToInteger

### ** Examples

library(oce)
twenty.five <- bcdToInteger(as.raw(0x25))
thirty.seven <- as.integer(as.raw(0x25))



cleanEx()
nameEx("beamUnspreadAdp")
### * beamUnspreadAdp

flush(stderr()); flush(stdout())

### Name: beamUnspreadAdp
### Title: Adjust ADP Signal for Spherical Spreading
### Aliases: beamUnspreadAdp

### ** Examples

library(oce)
data(adp)
plot(adp, which=5) # beam 1 echo intensity
adp.att <- beamUnspreadAdp(adp)
plot(adp.att, which=5) # beam 1 echo intensity
## Profiles
par(mar=c(4, 4, 1, 1))
a <- adp[["a", "numeric"]]             # second arg yields matrix return value
distance <- adp[["distance"]]
plot(apply(a,2,mean), distance, type='l', xlim=c(0,256))
lines(apply(a,2,median), distance, type='l',col='red')
legend("topright",lwd=1,col=c("black","red"),legend=c("original","attenuated"))
## Image
plot(adp.att, which="amplitude",col=oce.colorsViridis(100))




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("binApply1D")
### * binApply1D

flush(stderr()); flush(stdout())

### Name: binApply1D
### Title: Apply a function to vector data
### Aliases: binApply1D

### ** Examples

library(oce)
## salinity profile with median and quartile 1 and 3
data(ctd)
p <- ctd[["pressure"]]
S <- ctd[["salinity"]]
q1 <- binApply1D(p, S, pretty(p, 30), function(x) quantile(x, 1/4))
q3 <- binApply1D(p, S, pretty(p, 30), function(x) quantile(x, 3/4))
plotProfile(ctd, "salinity", col='gray', type='n')
polygon(c(q1$result, rev(q3$result)),
c(q1$xmids, rev(q1$xmids)), col='gray')
points(S, p, pch=20)




cleanEx()
nameEx("binApply2D")
### * binApply2D

flush(stderr()); flush(stdout())

### Name: binApply2D
### Title: Apply a function to matrix data
### Aliases: binApply2D

### ** Examples

library(oce)




cleanEx()
nameEx("binAverage")
### * binAverage

flush(stderr()); flush(stdout())

### Name: binAverage
### Title: Bin-average a vector y, based on x values
### Aliases: binAverage

### ** Examples

library(oce)
## A. fake linear data
x <- seq(0, 100, 1)
y <- 1 + 2 * x
plot(x, y, pch=1)
ba <- binAverage(x, y)
points(ba$x, ba$y, pch=3, col='red', cex=3)

## B. fake quadratic data
y <- 1 + x ^2
plot(x, y, pch=1)
ba <- binAverage(x, y)
points(ba$x, ba$y, pch=3, col='red', cex=3)

## C. natural data
data(co2)
plot(co2)
avg <- binAverage(time(co2), co2, 1950, 2000, 2)
points(avg$x, avg$y, col='red')




cleanEx()
nameEx("binMean1D")
### * binMean1D

flush(stderr()); flush(stdout())

### Name: binMean1D
### Title: Bin-average f=f(x)
### Aliases: binMean1D

### ** Examples

library(oce)
data(ctd)
z <- ctd[["z"]]
T <- ctd[["temperature"]]
plot(T, z)
TT <- binMean1D(z, T, seq(-100, 0, 1))
lines(TT$result, TT$xmids, col='red')




cleanEx()
nameEx("binMean2D")
### * binMean2D

flush(stderr()); flush(stdout())

### Name: binMean2D
### Title: Bin-average f=f(x,y)
### Aliases: binMean2D

### ** Examples

library(oce)
x <- runif(500)
y <- runif(500)
f <- x + y
xb <- seq(0, 1, 0.1)
yb <- seq(0, 1, 0.2)
m <- binMean2D(x, y, f, xb, yb)
plot(x, y)
contour(m$xmids, m$ymids, m$result, add=TRUE, levels=seq(0, 2, 0.5), labcex=1)




cleanEx()
nameEx("binmapAdp")
### * binmapAdp

flush(stderr()); flush(stdout())

### Name: binmapAdp
### Title: Bin-map an ADP object
### Aliases: binmapAdp

### ** Examples

## Not run: 
##D library(oce)
##D beam <- read.oce("/data/archive/sleiwex/2008/moorings/m09/adp/rdi_2615/raw/adp_rdi_2615.000",
##D                  from=as.POSIXct("2008-06-26", tz="UTC"),
##D                  to=as.POSIXct("2008-06-26 00:10:00", tz="UTC"),
##D                  longitude=-69.73433, latitude=47.88126)
##D beam2 <- binmapAdp(beam)
##D plot(enuToOther(toEnu(beam), heading=-31.5))
##D plot(enuToOther(toEnu(beam2), heading=-31.5))
##D plot(beam, which=5:8) # backscatter amplitude
##D plot(beam2, which=5:8)
## End(Not run)




cleanEx()
nameEx("byteToBinary")
### * byteToBinary

flush(stderr()); flush(stdout())

### Name: byteToBinary
### Title: Format bytes as binary defunct
### Aliases: byteToBinary

### ** Examples

library(oce)
## Note comparison with rawToBits():
a <- as.raw(0x0a)
byteToBinary(a, "big")        # "00001010"
as.integer(rev(rawToBits(a))) # 0 0 0 0 1 0 1 0



cleanEx()
nameEx("cm")
### * cm

flush(stderr()); flush(stdout())

### Name: cm
### Title: A Current Meter (cm) Object
### Aliases: cm

### ** Examples

library(oce)
data(cm)
summary(cm)
plot(cm)




cleanEx()
nameEx("coastlineCut")
### * coastlineCut

flush(stderr()); flush(stdout())

### Name: coastlineCut
### Title: Cut a Coastline Object at Specified Longitude
### Aliases: coastlineCut

### ** Examples





cleanEx()
nameEx("colormap")
### * colormap

flush(stderr()); flush(stdout())

### Name: colormap
### Title: Calculate color map
### Aliases: colormap

### ** Examples

library(oce)
## Example 1. color scheme for points on xy plot
x <- seq(0, 1, length.out=40)
y <- sin(2 * pi * x)
par(mar=c(3, 3, 1, 1))
mar <- par('mar') # prevent margin creep by drawPalette()
## First, default breaks
c <- colormap(y)
drawPalette(c$zlim, col=c$col, breaks=c$breaks)
plot(x, y, bg=c$zcol, pch=21, cex=1)
grid()
par(mar=mar)
## Second, 100 breaks, yielding a smoother palette
c <- colormap(y, breaks=100)
drawPalette(c$zlim, col=c$col, breaks=c$breaks)
plot(x, y, bg=c$zcol, pch=21, cex=1)
grid()
par(mar=mar)

## Not run: 
##D ## Example 2. topographic image with a standard color scheme
##D par(mfrow=c(1,1))
##D data(topoWorld)
##D cm <- colormap(name="gmt_globe")
##D imagep(topoWorld, breaks=cm$breaks, col=cm$col)
##D 
##D ## Example 3. topographic image with modified colors,
##D ## black for depths below 4km.
##D cm <- colormap(name="gmt_globe")
##D deep <- cm$x0 < -4000
##D cm$col0[deep] <- 'black'
##D cm$col1[deep] <- 'black'
##D cm <- colormap(x0=cm$x0, x1=cm$x1, col0=cm$col0, col1=cm$col1)
##D imagep(topoWorld, breaks=cm$breaks, col=cm$col)
##D 
##D ## Example 4. image of world topography with water colorized
##D ## smoothly from violet at 8km depth to blue
##D ## at 4km depth, then blending in 0.5km increments
##D ## to white at the coast, with tan for land.
##D cm <- colormap(x0=c(-8000, -4000,   0,  100),
##D                x1=c(-4000,     0, 100, 5000),
##D                col0=c("violet","blue","white","tan"),
##D                col1=c("blue","white","tan","yellow"))
##D lon <- topoWorld[['longitude']]
##D lat <- topoWorld[['latitude']]
##D z <- topoWorld[['z']]
##D imagep(lon, lat, z, breaks=cm$breaks, col=cm$col)
##D contour(lon, lat, z, levels=0, add=TRUE)
##D 
##D ## Example 5. visualize GMT style color map
##D cm <- colormap(name="gmt_globe", debug=4)
##D plot(seq_along(cm$x0), cm$x0, pch=21, bg=cm$col0)
##D grid()
##D points(seq_along(cm$x1), cm$x1, pch=21, bg=cm$col1)
##D 
##D ## Example 6. colfunction
##D cm <- colormap(c(0, 1))
##D x <- 1:10
##D y <- (x - 5.5)^2
##D z <- seq(0, 1, length.out=length(x))
##D drawPalette(colormap=cm)
##D plot(x, y, pch=21, bg=cm$colfunction(z), cex=3)
## End(Not run)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("concatenate-adp-method")
### * concatenate-adp-method

flush(stderr()); flush(stdout())

### Name: concatenate,adp-method
### Title: Concatenate adp objects
### Aliases: concatenate,adp-method

### ** Examples

## 1. Split, then recombine, a ctd object.
data(ctd)
ctd1 <- subset(ctd, scan <= median(ctd[["scan"]]))
ctd2 <- subset(ctd, scan > median(ctd[["scan"]]))
CTD <- concatenate(ctd1, ctd2)

## 2. Split, then recombine, an adp object.
data(adp)
midtime <- median(adp[["time"]])
adp1 <- subset(adp, time <= midtime)
adp2 <- subset(adp, time > midtime)
ADP <- concatenate(adp1, adp2)

## Not run: 
##D ## 3. Download two met files and combine them.
##D met1 <- read.met(download.met(id=6358, year=2003, month=8))
##D met2 <- read.met(download.met(id=6358, year=2003, month=9))
##D MET <- concatenate(met1, met2)
## End(Not run)




cleanEx()
nameEx("concatenate-oce-method")
### * concatenate-oce-method

flush(stderr()); flush(stdout())

### Name: concatenate,oce-method
### Title: Concatenate oce objects
### Aliases: concatenate,oce-method

### ** Examples

## 1. Split, then recombine, a ctd object.
data(ctd)
ctd1 <- subset(ctd, scan <= median(ctd[["scan"]]))
ctd2 <- subset(ctd, scan > median(ctd[["scan"]]))
CTD <- concatenate(ctd1, ctd2)

## 2. Split, then recombine, an adp object.
data(adp)
midtime <- median(adp[["time"]])
adp1 <- subset(adp, time <= midtime)
adp2 <- subset(adp, time > midtime)
ADP <- concatenate(adp1, adp2)

## Not run: 
##D ## 3. Download two met files and combine them.
##D met1 <- read.met(download.met(id=6358, year=2003, month=8))
##D met2 <- read.met(download.met(id=6358, year=2003, month=9))
##D MET <- concatenate(met1, met2)
## End(Not run)




cleanEx()
nameEx("coriolis")
### * coriolis

flush(stderr()); flush(stdout())

### Name: coriolis
### Title: Coriolis parameter on rotating earth
### Aliases: coriolis

### ** Examples

C <- coriolis(45) # 1e-4



cleanEx()
nameEx("ctd-class")
### * ctd-class

flush(stderr()); flush(stdout())

### Name: ctd-class
### Title: Class to Store CTD (or general hydrographic) Data
### Aliases: ctd-class

### ** Examples


# 1. Create a ctd object with fake data.
a <- as.ctd(salinity=35+1:3/10, temperature=10-1:3/10, pressure=1:3)
summary(a)

# 2. Fix a typo in a station latitude (fake! it's actually okay)
data(ctd)
ctd <- oceSetMetadata(ctd, "latitude", ctd[["latitude"]]-0.001,
                     "fix latitude typo in log book")




cleanEx()
nameEx("ctd")
### * ctd

flush(stderr()); flush(stdout())

### Name: ctd
### Title: A CTD profile in Halifax Harbour
### Aliases: ctd

### ** Examples

## Not run: 
##D library(oce)
##D data(ctd)
##D plot(ctd)
## End(Not run)




cleanEx()
nameEx("ctd.cnv")
### * ctd.cnv

flush(stderr()); flush(stdout())

### Name: ctd.cnv
### Title: Sample ctd dataset in .cnv format
### Aliases: ctd.cnv

### ** Examples

## Not run: 
##D read.oce(system.file("extdata", "ctd.cnv", package="oce"))
## End(Not run)




cleanEx()
nameEx("ctdDecimate")
### * ctdDecimate

flush(stderr()); flush(stdout())

### Name: ctdDecimate
### Title: Decimate a CTD profile
### Aliases: ctdDecimate

### ** Examples

library(oce)
data(ctd)
plotProfile(ctd, "salinity", ylim=c(10, 0))
p <- seq(0, 45, 1)
ctd2 <- ctdDecimate(ctd, p=p)
lines(ctd2[["salinity"]], ctd2[["pressure"]], col="blue")
p <- seq(0, 45, 1)
ctd3 <- ctdDecimate(ctd, p=p, method=function(x, y, xout)
                    predict(smooth.spline(x, y, df=30), xout)$y)
lines(ctd3[["salinity"]], ctd3[["pressure"]], col="red")





cleanEx()
nameEx("ctdFindProfiles")
### * ctdFindProfiles

flush(stderr()); flush(stdout())

### Name: ctdFindProfiles
### Title: Find Profiles within a Tow-Yow CTD Record
### Aliases: ctdFindProfiles

### ** Examples

library(oce)

## Not run: 
##D # Example 1. Find profiles within a towyo file, as can result
##D # if the CTD is cycled within the water column as the ship
##D # moves.
##D profiles <- ctdFindProfiles(towyo)
## End(Not run)

## Not run: 
##D # Example 2. Use a moving average to smooth pressure, instead of the
##D # default smooth.spline() method. This might avoid a tendency of
##D # the default scheme to miss some profiles in a long towyo.
##D movingAverage <- function(x, n = 11, ...)
##D {
##D     f <- rep(1/n, n)
##D     stats::filter(x, f, ...)
##D }
##D casts <- ctdFindProfiles(towyo, smoother=movingAverage)
## End(Not run)

## Not run: 
##D # Example 3: glider data read into a ctd object. Chop
##D # into profiles by looking for pressure jumps exceeding
##D # 10 dbar.
##D breaks <- which(diff(gliderAsCtd[["pressure"]]) > 10)
##D profiles <- ctdFindProfiles(gliderAsCtd, breaks=breaks)
## End(Not run)




cleanEx()
nameEx("ctdTrim")
### * ctdTrim

flush(stderr()); flush(stdout())

### Name: ctdTrim
### Title: Trim Beginning and Ending of a CTD cast
### Aliases: ctdTrim

### ** Examples

## Not run: 
##D library(oce)
##D data(ctdRaw)
##D plot(ctdRaw) # barely recognizable, due to pre- and post-cast junk
##D plot(ctdTrim(ctdRaw)) # looks like a real profile ...
##D plot(ctdDecimate(ctdTrim(ctdRaw),method="boxcar")) # ... smoothed
##D # Demonstrate use of a function. The scan limits were chosen
##D # by using locator(2) on a graph made by plotScan(ctdRaw).
##D trimByIndex<-function(data, parameters) {
##D   parameters[1] < data$scan & data$scan < parameters[2]
##D }
##D trimmed <- ctdTrim(ctdRaw, trimByIndex, parameters=c(130, 380))
##D plot(trimmed)
## End(Not run)




cleanEx()
nameEx("ctimeToSeconds")
### * ctimeToSeconds

flush(stderr()); flush(stdout())

### Name: ctimeToSeconds
### Title: Interpret a character string as a time interval
### Aliases: ctimeToSeconds

### ** Examples

library(oce)
cat("10      = ", ctimeToSeconds("10"), "s\n", sep="")
cat("01:04   = ", ctimeToSeconds("01:04"), "s\n", sep="")
cat("1:00:00 = ", ctimeToSeconds("1:00:00"), "s\n", sep="")



cleanEx()
nameEx("curl")
### * curl

flush(stderr()); flush(stdout())

### Name: curl
### Title: Curl of 2D vector field
### Aliases: curl

### ** Examples

library(oce)
## 1. Shear flow with uniform curl.
x <- 1:4
y <- 1:10
u <- outer(x, y, function(x, y) y/2)
v <- outer(x, y, function(x, y) -x/2)
C <- curl(u, v, x, y, FALSE)

## 2. Rankine vortex: constant curl inside circle, zero outside
rankine <- function(x, y)
{
    r <- sqrt(x^2 + y^2)
    theta <- atan2(y, x)
    speed <- ifelse(r < 1, 0.5*r, 0.5/r)
    list(u=-speed*sin(theta), v=speed*cos(theta))
}
x <- seq(-2, 2, length.out=100)
y <- seq(-2, 2, length.out=50)
u <- outer(x, y, function(x, y) rankine(x, y)$u)
v <- outer(x, y, function(x, y) rankine(x, y)$v)
C <- curl(u, v, x, y, FALSE)
## plot results
par(mfrow=c(2, 2))
imagep(x, y, u, zlab="u", asp=1)
imagep(x, y, v, zlab="v", asp=1)
imagep(x, y, C$curl, zlab="curl", asp=1)
hist(C$curl, breaks=100)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("d200321-001.ctd")
### * d200321-001.ctd

flush(stderr()); flush(stdout())

### Name: d200321-001.ctd
### Title: Sample ctd dataset in .ctd format
### Aliases: d200321-001.ctd

### ** Examples

## Not run: 
##D read.oce(system.file("extdata", "d200321-001.ctd", package="oce"))
## End(Not run)




cleanEx()
nameEx("d201211_0011.cnv")
### * d201211_0011.cnv

flush(stderr()); flush(stdout())

### Name: d201211_0011.cnv
### Title: Sample ctd dataset in .cnv format
### Aliases: d201211_0011.cnv

### ** Examples

## Not run: 
##D read.oce(system.file("extdata", "d201211_0011.cnv", package="oce"))
## End(Not run)




cleanEx()
nameEx("dataLabel")
### * dataLabel

flush(stderr()); flush(stdout())

### Name: dataLabel
### Title: Try to associate data names with units, for use by summary()
### Aliases: dataLabel

### ** Examples

library(oce)
data(ctd)
dataLabel(names(ctd@data), ctd@metadata$units)



cleanEx()
nameEx("decimate")
### * decimate

flush(stderr()); flush(stdout())

### Name: decimate
### Title: Smooth and Decimate, or Subsample, an Oce Object
### Aliases: decimate

### ** Examples

library(oce)
data(adp)
plot(adp)
adpDec <- decimate(adp,by=2,filter=c(1/4, 1/2, 1/4))
plot(adpDec)



cleanEx()
nameEx("decodeTime")
### * decodeTime

flush(stderr()); flush(stdout())

### Name: decodeTime
### Title: Oce Version of as.POSIXct
### Aliases: decodeTime

### ** Examples

decodeTime("July 1 2013 01:02:03")
decodeTime("Jul 1 2013 01:02:03")
decodeTime("1 July 2013 01:02:03")
decodeTime("1 Jul 2013 01:02:03")
decodeTime("2013-07-01 01:02:03")
decodeTime("2013/07/01 01:02:03")
decodeTime("2013/07/01")




cleanEx()
nameEx("despike")
### * despike

flush(stderr()); flush(stdout())

### Name: despike
### Title: Remove spikes from a time series
### Aliases: despike

### ** Examples

n <- 50
x <- 1:n
y <- rnorm(n=n)
y[n/2] <- 10                    # 10 standard deviations
plot(x, y, type='l')
lines(x, despike(y), col='red')
lines(x, despike(y, reference="smooth"), col='darkgreen')
lines(x, despike(y, reference="trim", min=-3, max=3), col='blue')
legend("topright", lwd=1, col=c("black", "red", "darkgreen", "blue"),
       legend=c("raw", "median", "smooth", "trim"))

# add a spike to a CTD object
data(ctd)
plot(ctd)
T <- ctd[["temperature"]]
T[10] <- T[10] + 10
ctd[["temperature"]] <- T
CTD <- despike(ctd)
plot(CTD)



cleanEx()
nameEx("detrend")
### * detrend

flush(stderr()); flush(stdout())

### Name: detrend
### Title: Detrend a set of observations
### Aliases: detrend

### ** Examples

x <- seq(0, 0.9 * pi, length.out=50)
y <- sin(x)
y[1] <- NA
y[10] <- NA
plot(x, y, ylim=c(0, 1))
d <- detrend(x, y)
points(x, d$Y, pch=20)
abline(d$a, d$b, col='blue')
abline(h=0)
points(x, d$Y + d$a + d$b * x, col='blue', pch='+')



cleanEx()
nameEx("download.amsr")
### * download.amsr

flush(stderr()); flush(stdout())

### Name: download.amsr
### Title: Download and Cache an amsr File
### Aliases: download.amsr

### ** Examples

## Not run: 
##D ## The download takes several seconds.
##D f <- download.amsr(2017, 1, 14) # Jan 14, 2017
##D d <- read.amsr(f)
##D plot(d)
##D mtext(d[["filename"]], side=3, line=0, adj=0)
## End(Not run)




cleanEx()
nameEx("download.coastline")
### * download.coastline

flush(stderr()); flush(stdout())

### Name: download.coastline
### Title: Download a coastline File
### Aliases: download.coastline

### ** Examples

## Not run: 
##D library(oce)
##D # User must create directory ~/data/coastline first.
##D # As of September 2016, the downloaded file, named
##D # "ne_50m_coastline.zip", occupies 443K bytes.
##D filename <- download.coastline(destdir="~/data/coastline")
##D coastline <- read.coastline(filename)
##D plot(coastline)
## End(Not run)




cleanEx()
nameEx("download.met")
### * download.met

flush(stderr()); flush(stdout())

### Name: download.met
### Title: Download and Cache a met File
### Aliases: download.met

### ** Examples

## Not run: 
##D library(oce)
##D ## Download data for Halifax International Airport, in September
##D ## of 2003. (This dataset is used for data(met) provided with oce.)
##D metFile <- download.met(6358, 2003, 9, destdir=".")
##D met <- read.met(metFile)
## End(Not run)




cleanEx()
nameEx("download.topo")
### * download.topo

flush(stderr()); flush(stdout())

### Name: download.topo
### Title: Download and Cache a topo File
### Aliases: download.topo

### ** Examples

## Not run: 
##D library(oce)
##D topoFile <- download.topo(west=-66, east=-60, south=43, north=47,
##D                           resolution=1, destdir="~/data/topo")
##D topo <- read.topo(topoFile)
##D imagep(topo, zlim=c(-400, 400), col=oceColorsTwo, drawTriangles=TRUE)
##D if (requireNamespace("ocedata", quietly=TRUE)) {
##D     data(coastlineWorldFine, package="ocedata")
##D     lines(coastlineWorldFine[["longitude"]], coastlineWorldFine[["latitude"]])
##D }
## End(Not run)




cleanEx()
nameEx("drawDirectionField")
### * drawDirectionField

flush(stderr()); flush(stdout())

### Name: drawDirectionField
### Title: Draw a Direction Field
### Aliases: drawDirectionField

### ** Examples

library(oce)
plot(c(-1.5, 1.5), c(-1.5, 1.5), xlab="", ylab="", type='n')
drawDirectionField(x=rep(0, 2), y=rep(0, 2), u=c(1, 1), v=c(1, -1), scalex=0.5, add=TRUE)
plot(c(-1.5, 1.5), c(-1.5, 1.5), xlab="", ylab="", type='n')
drawDirectionField(x=rep(0, 2), y=rep(0, 2), u=c(1, 1), v=c(1, -1), scalex=0.5, add=TRUE,
                   type=2)

## 2D example
x <- seq(-2, 2, 0.1)
y <- x
xx <- expand.grid(x, y)[,1]
yy <- expand.grid(x, y)[,2]
z <- matrix(xx*exp(-xx^2 -yy^2), nrow=length(x))
gz <- grad(z, x, y)
drawDirectionField(x, y, gz$gx, gz$gy, scalex=0.5, type=2, len=0.02)
oceContour(x, y, z, add=TRUE)




cleanEx()
nameEx("drawPalette")
### * drawPalette

flush(stderr()); flush(stdout())

### Name: drawPalette
### Title: Draw a palette, leaving margins suitable for accompanying plot
### Aliases: drawPalette

### ** Examples


library(oce)
par(mgp=getOption("oceMgp"))

## 1. A three-panel plot
par(mfrow=c(3, 1), mar=c(3, 3, 1, 1))
omar <- par('mar')                 # save initial margin

## 1a. top panel: simple case with Viridis scheme
drawPalette(zlim=c(0, 1), col=oce.colorsViridis(10))
plot(1:10, 1:10, col=oce.colorsViridis(10)[1:10],pch=20,cex=3,xlab='x', ylab='y')
par(mar=omar)                      # reset margin

## 1b. middle panel: colormap
cm <- colormap(name="gmt_globe")
drawPalette(colormap=cm)
icol <- seq_along(cm$col)
plot(icol, cm$breaks[icol], pch=20, cex=2, col=cm$col,
     xlab="Palette index", ylab="Palette breaks")
par(mar=omar)                      # reset margin

## 1c. bottom panel: space for palette (to line up graphs)
drawPalette(plot=FALSE)
plot(1:10, 1:10, col=oce.colorsViridis(10)[1:10],pch=20,cex=3,xlab='x', ylab='y')
par(mar=omar)                      # reset margin

# 2. Use layout to mimic the action of imagep(), with the width
# of the palette region being 14 percent of figure width.
d <- 0.14
layout(matrix(1:2,nrow=1), widths=c(1-d,d))
image(volcano, col=oce.colorsViridis(100), zlim=c(90, 200))
contour(volcano, add=TRUE)
drawPalette(c(90, 200), fullpage=TRUE, col=oce.colorsViridis)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("enuToOtherAdp")
### * enuToOtherAdp

flush(stderr()); flush(stdout())

### Name: enuToOtherAdp
### Title: Convert ADP ENU to Rotated Coordinate
### Aliases: enuToOtherAdp

### ** Examples


library(oce)
data(adp)
o <- enuToOtherAdp(adp, heading=-31.5)
plot(o, which=1:3)




cleanEx()
nameEx("errorbars")
### * errorbars

flush(stderr()); flush(stdout())

### Name: errorbars
### Title: Draw error bars on an existing xy diagram
### Aliases: errorbars

### ** Examples

library(oce)
data(ctd)
S <- ctd[["salinity"]]
T <- ctd[["temperature"]]
plot(S, T)
errorbars(S, T, 0.05, 0.5)



cleanEx()
nameEx("fillGap")
### * fillGap

flush(stderr()); flush(stdout())

### Name: fillGap
### Title: Fill a gap in an oce object
### Aliases: fillGap

### ** Examples

library(oce)
# Integers
x <- c(1:2, NA, NA, 5:6)
y <- fillGap(x)
print(data.frame(x,y))
# Floats
x <- x + 0.1
y <- fillGap(x)
print(data.frame(x,y))



cleanEx()
nameEx("formatCI")
### * formatCI

flush(stderr()); flush(stdout())

### Name: formatCI
### Title: Confidence interval in parenthetic notation
### Aliases: formatCI

### ** Examples

x <- seq(0, 1, length.out=300)
y <- rnorm(n=300, mean=10, sd=1) * x
m <- lm(y~x)
print(formatCI(model=m))



cleanEx()
nameEx("formatPosition")
### * formatPosition

flush(stderr()); flush(stdout())

### Name: formatPosition
### Title: Format Geographical Position in Degrees and Minutes
### Aliases: formatPosition

### ** Examples

library(oce)
formatPosition(10+1:10/60+2.8/3600)
formatPosition(10+1:10/60+2.8/3600, type="string")




cleanEx()
nameEx("geodDist")
### * geodDist

flush(stderr()); flush(stdout())

### Name: geodDist
### Title: Compute Geodesic Distance on Surface of Earth
### Aliases: geodDist

### ** Examples

library(oce)
km <- geodDist(100, 45, 100, 46)
data(section)
geodDist(section)
geodDist(section, alongPath=TRUE)




cleanEx()
nameEx("geodGc")
### * geodGc

flush(stderr()); flush(stdout())

### Name: geodGc
### Title: Great-circle Segments Between Points on Earth
### Aliases: geodGc

### ** Examples





cleanEx()
nameEx("geodXy")
### * geodXy

flush(stderr()); flush(stdout())

### Name: geodXy
### Title: Convert From Geographical to Geodesic Coordinates
### Aliases: geodXy

### ** Examples




cleanEx()
nameEx("grad")
### * grad

flush(stderr()); flush(stdout())

### Name: grad
### Title: Calculate Matrix Gradient
### Aliases: grad

### ** Examples

## 1. Built-in volcano dataset
g <- grad(volcano)
par(mfrow=c(2, 2), mar=c(3, 3, 1, 1), mgp=c(2, 0.7, 0))
imagep(volcano, zlab="h")
imagep(g$g, zlab="|grad(h)|")
zlim <- c(-1, 1) * max(g$g)
imagep(g$gx, zlab="dh/dx", zlim=zlim)
imagep(g$gy, zlab="dh/dy", zlim=zlim)

## 2. Geostrophic flow around an eddy
library(oce)
dx <- 5e3
dy <- 10e3
x <- seq(-200e3, 200e3, dx)
y <- seq(-200e3, 200e3, dy)
R <- 100e3
h <- outer(x, y, function(x, y) 500*exp(-(x^2+y^2)/R^2))
grad <- grad(h, x, y)
par(mfrow=c(2, 2), mar=c(3, 3, 1, 1), mgp=c(2, 0.7, 0))
contour(x,y,h,asp=1, main=expression(h))
f <- 1e-4
gprime <- 9.8 * 1 / 1024
u <- -(gprime / f) * grad$gy
v <-  (gprime / f) * grad$gx
contour(x, y, u, asp=1, main=expression(u))
contour(x, y, v, asp=1, main=expression(v))
contour(x, y, sqrt(u^2+v^2), asp=1, main=expression(speed))




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("gravity")
### * gravity

flush(stderr()); flush(stdout())

### Name: gravity
### Title: Acceleration due to earth gravity
### Aliases: gravity

### ** Examples

g <- gravity(45) # 9.8



cleanEx()
nameEx("handleFlags-adp-method")
### * handleFlags-adp-method

flush(stderr()); flush(stdout())

### Name: handleFlags,adp-method
### Title: Handle Flags in adp Objects
### Aliases: handleFlags,adp-method

### ** Examples

# Flag low "goodness" or high "error beam" values.
library(oce)
data(adp)
# Same as Example 2 of ?'setFlags,adp-method'
v <- adp[["v"]]
i2 <- array(FALSE, dim=dim(v))
g <- adp[["g", "numeric"]]
# Thresholds on percent "goodness" and error "velocity"
G <- 25
V4 <- 0.45
for (k in 1:3)
    i2[,,k] <- ((g[,,k]+g[,,4]) < G) | (v[,,4] > V4)
adpQC <- initializeFlags(adp, "v", 2)
adpQC <- setFlags(adpQC, "v", i2, 3)
adpClean <- handleFlags(adpQC, flags=list(3), actions=list("NA"))
# Demonstrate (subtle) change graphically.
par(mfcol=c(2, 1))
plot(adp, which="u1")
plot(adpClean, which="u1")




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("handleFlags-argo-method")
### * handleFlags-argo-method

flush(stderr()); flush(stdout())

### Name: handleFlags,argo-method
### Title: Handle Flags in ARGO Objects
### Aliases: handleFlags,argo-method handleFlags.argo

### ** Examples

library(oce)
data(argo)
argoNew <- handleFlags(argo)
# Demonstrate replacement, looking at the second profile
f <- argo[["salinityFlag"]][,2]
df <- data.frame(flag=f, orig=argo[["salinity"]][,2], new=argoNew[["salinity"]][,2])
df[11:15,] # notice line 13




cleanEx()
nameEx("handleFlags-ctd-method")
### * handleFlags-ctd-method

flush(stderr()); flush(stdout())

### Name: handleFlags,ctd-method
### Title: Handle Flags in ctd Objects
### Aliases: handleFlags,ctd-method

### ** Examples

library(oce)
data(section)
stn <- section[["station", 100]]
# 1. Default: anything not flagged as 2 is set to NA, to focus
# solely on 'good', in the World Hydrographic Program scheme.
STN1 <- handleFlags(stn, flags=list(c(1, 3:9)))
data.frame(old=stn[["salinity"]], new=STN1[["salinity"]], salinityFlag=stn[["salinityFlag"]])

# 2. Use bottle salinity, if it is good and ctd is bad
replace <- 2 == stn[["salinityBottleFlag"]] && 2 != stn[["salinityFlag"]]
S <- ifelse(replace, stn[["salinityBottle"]], stn[["salinity"]])
STN2 <- oceSetData(stn, "salinity", S)

# 3. Use smoothed TS relationship to nudge questionable data.
f <- function(x) {
  S <- x[["salinity"]]
  T <- x[["temperature"]]
  df <- 0.5 * length(S) # smooths a bit
  sp <- smooth.spline(T, S, df=df)
  0.5 * (S + predict(sp, T)$y)
}
par(mfrow=c(1,2))
STN3 <- handleFlags(stn, flags=list(salinity=c(1,3:9)), action=list(salinity=f))
plotProfile(stn, "salinity", mar=c(3, 3, 3, 1))
p <- stn[["pressure"]]
par(mar=c(3, 3, 3, 1))
plot(STN3[["salinity"]] - stn[["salinity"]], p, ylim=rev(range(p)))

# 4. Single-variable flags (vector specification)
data(section)
# Multiple-flag scheme: one per data item
A <- section[["station", 100]]
deep <- A[["pressure"]] > 1500
flag <- ifelse(deep, 7, 2)
for (flagName in names(A[["flags"]]))
    A[[paste(flagName, "Flag", sep="")]] <- flag
Af <- handleFlags(A)
expect_equal(is.na(Af[["salinity"]]), deep)

# 5. Single-variable flags (list specification)
B <- section[["station", 100]]
B[["flags"]] <- list(flag)
Bf <- handleFlags(B)
expect_equal(is.na(Bf[["salinity"]]), deep)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("handleFlags-section-method")
### * handleFlags-section-method

flush(stderr()); flush(stdout())

### Name: handleFlags,section-method
### Title: Handle flags in Section Objects
### Aliases: handleFlags,section-method handleFlags.section

### ** Examples

library(oce)
data(section)
section2 <- handleFlags(section, flags=c(1,3:9))
par(mfrow=c(2, 1))
plotTS(section)
plotTS(section2)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("imagep")
### * imagep

flush(stderr()); flush(stdout())

### Name: imagep
### Title: Plot an Image with a Color Palette
### Aliases: imagep

### ** Examples

library(oce)

# 1. simplest use
imagep(volcano)

# 2. something oceanographic (internal-wave speed)
h <- seq(0, 50, length.out=100)
drho <- seq(1, 3, length.out=200)
speed <- outer(h, drho, function(drho, h) sqrt(9.8 * drho * h / 1024))
imagep(h, drho, speed, xlab="Equivalent depth [m]",
ylab=expression(paste(Delta*rho, " [kg/m^3]")),
zlab="Internal-wave speed [m/s]")

# 3. fancy labelling on atan() function
x <- seq(0, 1, 0.01)
y <- seq(0, 1, 0.01)
angle <- outer(x,y,function(x,y) atan2(y,x))
imagep(x, y, angle, filledContour=TRUE, breaks=c(0, pi/4, pi/2),
       col=c("lightgray", "darkgray"),
       at=c(0, pi/4, pi/2),
       labels=c(0, expression(pi/4), expression(pi/2)))

# 5. y-axis flipping
par(mfrow=c(2,2))
data(adp)
d <- adp[["distance"]]
t <- adp[["time"]]
u <- adp[["v"]][ , ,1]
imagep(t, d, u, drawTimeRange=FALSE)
mtext("normal")
imagep(t, d, u, flipy=TRUE, drawTimeRange=FALSE)
mtext("flipy")
imagep(t, d, u, ylim=rev(range(d)), drawTimeRange=FALSE)
mtext("ylim")
imagep(t, d, u, ylim=rev(range(d)), flipy=TRUE, drawTimeRange=FALSE)
mtext("flipy and ylim")
par(mfrow=c(1,1))

# 6. a colormap case
data(topoWorld)
cm <- colormap(name="gmt_globe")
imagep(topoWorld, colormap=cm)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("initialize-ctd-method")
### * initialize-ctd-method

flush(stderr()); flush(stdout())

### Name: initialize,ctd-method
### Title: Initialize storage for a ctd object
### Aliases: initialize,ctd-method

### ** Examples


## 1. empty
new("ctd")

## 2. fake data with no location information, so can only
##    plot with the UNESCO equation of state.
##    NOTE: always name arguments, in case the default order gets changed
ctd <- new("ctd", salinity=35+1:3/10, temperature=10-1:3/10, pressure=1:3)
summary(ctd)
plot(ctd, eos="unesco")

## 3. as 2, but insert location and plot with GSW equation of state.
ctd <- oceSetMetadata(ctd, "latitude", 44)
ctd <- oceSetMetadata(ctd, "longitude", -63)
plot(ctd, eos="gsw")




cleanEx()
nameEx("initializeFlagScheme-section-method")
### * initializeFlagScheme-section-method

flush(stderr()); flush(stdout())

### Name: initializeFlagScheme,section-method
### Title: Establish a data-quality scheme for a section object
### Aliases: initializeFlagScheme,section-method

### ** Examples

## Not run: 
##D data(section)
##D section <- read.section("a03_hy1.csv", sectionId="a03", institute="SIO",
##D                         ship="R/V Professor Multanovskiy", scientist="Vladimir Tereschenkov")
##D sectionWithFlags <- initializeFlagScheme(section, "WHP bottle")
##D station1 <- sectionWithFlags[["station", 1]]
##D str(station1[["flagScheme"]])
## End(Not run)



cleanEx()
nameEx("integerToAscii")
### * integerToAscii

flush(stderr()); flush(stdout())

### Name: integerToAscii
### Title: Decode integer to corresponding ASCII code
### Aliases: integerToAscii

### ** Examples

library(oce)
A <- integerToAscii(65)
cat("A=", A, "\n")



cleanEx()
nameEx("integrateTrapezoid")
### * integrateTrapezoid

flush(stderr()); flush(stdout())

### Name: integrateTrapezoid
### Title: Trapezoidal Integration
### Aliases: integrateTrapezoid

### ** Examples

x <- seq(0, 1, length.out=10) # try larger length.out to see if area approaches 2
y <- 2*x + 3*x^2
A <- integrateTrapezoid(x, y)
dA <- integrateTrapezoid(x, y, "dA")
cA <- integrateTrapezoid(x, y, "cA")
print(A)
print(sum(dA))
print(tail(cA, 1))
print(integrateTrapezoid(diff(x[1:2]), y))
print(integrateTrapezoid(y))



cleanEx()
nameEx("interpBarnes")
### * interpBarnes

flush(stderr()); flush(stdout())

### Name: interpBarnes
### Title: Grid data using Barnes algorithm
### Aliases: interpBarnes

### ** Examples

library(oce)

# 1. contouring example, with wind-speed data from Koch et al. (1983)
data(wind)
u <- interpBarnes(wind$x, wind$y, wind$z)
contour(u$xg, u$yg, u$zg, labcex=1)
text(wind$x, wind$y, wind$z, cex=0.7, col="blue")
title("Numbers are the data")

# 2. As 1, but blank out spots where data are sparse
u <- interpBarnes(wind$x, wind$y, wind$z, trim=0.1)
contour(u$xg, u$yg, u$zg, level=seq(0, 30, 1))
points(wind$x, wind$y, cex=1.5, pch=20, col="blue")

# 3. As 1, but interpolate back to points, and display the percent mismatch
u <- interpBarnes(wind$x, wind$y, wind$z)
contour(u$xg, u$yg, u$zg, labcex=1)
mismatch <- 100 * (wind$z - u$zd) / wind$z
text(wind$x, wind$y, round(mismatch), col="blue")
title("Numbers are percent mismatch between grid and data")


# 4. As 3, but contour the mismatch
mismatchGrid <- interpBarnes(wind$x, wind$y, mismatch)
contour(mismatchGrid$xg, mismatchGrid$yg, mismatchGrid$zg, labcex=1)

# 5. One-dimensional example, smoothing a salinity profile
data(ctd)
p <- ctd[["pressure"]]
y <- rep(1, length(p)) # fake y data, with arbitrary value
S <- ctd[["salinity"]]
pg <- pretty(p, n=100)
g <- interpBarnes(p, y, S, xg=pg, xr=1)
plot(S, p, cex=0.5, col="blue", ylim=rev(range(p)))
lines(g$zg, g$xg, col="red")



cleanEx()
nameEx("julianCenturyAnomaly")
### * julianCenturyAnomaly

flush(stderr()); flush(stdout())

### Name: julianCenturyAnomaly
### Title: Julian-Day number to Julian century
### Aliases: julianCenturyAnomaly

### ** Examples


t <- ISOdatetime(1978, 11, 13, 4, 35, 0, tz="UTC")
jca <- julianCenturyAnomaly(julianDay(t))
cat(format(t), "is Julian Century anomaly", format(jca, digits=8), "\n")




cleanEx()
nameEx("julianDay")
### * julianDay

flush(stderr()); flush(stdout())

### Name: julianDay
### Title: Convert a time to a Julian day
### Aliases: julianDay

### ** Examples

## example from Meeus
t <- ISOdatetime(1977, 4, 26, hour=0, min=0, sec=0, tz="UTC")+0.4*86400
expect_equal(julianDay(t), 2443259.9)




cleanEx()
nameEx("lobo")
### * lobo

flush(stderr()); flush(stdout())

### Name: lobo
### Title: Sample LOBO Dataset
### Aliases: lobo

### ** Examples

library(oce)
data(lobo)
summary(lobo)
plot(lobo)




cleanEx()
nameEx("lon360")
### * lon360

flush(stderr()); flush(stdout())

### Name: lon360
### Title: Alter longitudes from -180:180 to 0:360 convention
### Aliases: lon360

### ** Examples

lon360(c(179, -179))



cleanEx()
nameEx("lonlat2map")
### * lonlat2map

flush(stderr()); flush(stdout())

### Name: lonlat2map
### Title: Convert Longitude and Latitude to X and Y
### Aliases: lonlat2map

### ** Examples

library(oce)
## Cape Split, in the Minas Basin of the Bay of Fundy
cs <- list(longitude=-64.49657, latitude=45.33462)
xy <- lonlat2map(cs, projection="+proj=merc")
map2lonlat(xy)




cleanEx()
nameEx("lonlat2utm")
### * lonlat2utm

flush(stderr()); flush(stdout())

### Name: lonlat2utm
### Title: Convert Longitude and Latitude to UTM
### Aliases: lonlat2utm

### ** Examples

library(oce)
## Cape Split, in the Minas Basin of the Bay of Fundy
lonlat2utm(-64.496567, 45.334626)




cleanEx()
nameEx("lookWithin")
### * lookWithin

flush(stderr()); flush(stdout())

### Name: lookWithin
### Title: Look Within the First Element of a List for Replacement Values
### Aliases: lookWithin

### ** Examples

## 1. If first item is not a CTD object, just return the input
lookWithin(list(a=1, b=2)) # returns a list
## 2. Extract salinity from a CTD object
data(ctd)
str(lookWithin(list(salinity=ctd)))
## 3. Extract salinity and temperature. Note that the
## value specified for temperature is ignored; all that matters
## is that temperature is named.
str(lookWithin(list(salinity=ctd, temperature=NULL)))
## 4. How it is used by swRho()
rho1 <- swRho(ctd, eos="unesco")
rho2 <- swRho(ctd[["salinity"]], ctd[["temperature"]], ctd[["pressure"]], eos="unesco")
expect_equal(rho1, rho2)



cleanEx()
nameEx("lowpass")
### * lowpass

flush(stderr()); flush(stdout())

### Name: lowpass
### Title: Perform lowpass digital filtering
### Aliases: lowpass

### ** Examples

library(oce)
par(mfrow=c(1, 2), mar=c(4, 4, 1, 1))
coef <- lowpass(n=5, coefficients=TRUE)
plot(-2:2, coef, ylim=c(0, 1), xlab="Lag", ylab="Coefficient")
x <- seq(-5, 5) + rnorm(11)
plot(1:11, x, type='o', xlab="time", ylab="x and X")
X <- lowpass(x, n=5)
lines(1:11, X, col=2)
points(1:11, X, col=2)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("magneticField")
### * magneticField

flush(stderr()); flush(stdout())

### Name: magneticField
### Title: Earth magnetic declination, inclination, and intensity
### Aliases: magneticField

### ** Examples

library(oce)
# 1. Today's value at Halifax NS
magneticField(-(63+36/60), 44+39/60, Sys.Date())

# 2. World map of declination in year 2000.

# 3. Declination differences between versions 12 and 13



cleanEx()
nameEx("makeFilter")
### * makeFilter

flush(stderr()); flush(stdout())

### Name: makeFilter
### Title: Make a digital filter
### Aliases: makeFilter

### ** Examples

library(oce)

# 1. Demonstrate step-function response
y <- c(rep(1, 10), rep(-1, 10))
x <- seq_along(y)
plot(x, y, type='o', ylim=c(-1.05, 1.05))
BH <- makeFilter("blackman-harris", 11, asKernel=FALSE)
H <- makeFilter("hamming", 11, asKernel=FALSE)
yBH <- stats::filter(y, BH)
points(x, yBH, col=2, type='o')
yH <- stats::filter(y, H)
points(yH, col=3, type='o')
legend("topright", col=1:3, cex=2/3, pch=1,
       legend=c("input", "Blackman Harris", "Hamming"))

# 2. Show theoretical and practical filter gain, where
#    the latter is based on random white noise, and
#    includes a particular value for the spans
#    argument of spectrum(), etc.
## Not run: 
##D # need signal package for this example
##D r <- rnorm(2048)
##D rh <- stats::filter(r, H)
##D rh <- rh[is.finite(rh)] # kludge to remove NA at start/end
##D sR <- spectrum(r, plot=FALSE, spans=c(11, 5, 3))
##D sRH <- spectrum(rh, plot=FALSE, spans=c(11, 5, 3))
##D par(mfrow=c(2, 1), mar=c(3, 3, 1, 1), mgp=c(2, 0.7, 0))
##D plot(sR$freq, sRH$spec/sR$spec, xlab="Frequency", ylab="Power Transfer",
##D      type='l', lwd=5, col='gray')
##D theory <- freqz(H, n=seq(0,pi,length.out=100))
##D # Note we must square the modulus for the power spectrum
##D lines(theory$f/pi/2, Mod(theory$h)^2, lwd=1, col='red')
##D grid()
##D legend("topright", col=c("gray", "red"), lwd=c(5, 1), cex=2/3,
##D        legend=c("Practical", "Theory"), bg="white")
##D plot(log10(sR$freq), log10(sRH$spec/sR$spec),
##D      xlab="log10 Frequency", ylab="log10 Power Transfer",
##D      type='l', lwd=5, col='gray')
##D theory <- freqz(H, n=seq(0,pi,length.out=100))
##D # Note we must square the modulus for the power spectrum
##D lines(log10(theory$f/pi/2), log10(Mod(theory$h)^2), lwd=1, col='red')
##D grid()
##D legend("topright", col=c("gray", "red"), lwd=c(5, 1), cex=2/3,
##D        legend=c("Practical", "Theory"), bg="white")
## End(Not run)



cleanEx()
nameEx("map2lonlat")
### * map2lonlat

flush(stderr()); flush(stdout())

### Name: map2lonlat
### Title: Convert X and Y to Longitude and Latitude
### Aliases: map2lonlat

### ** Examples

library(oce)
## Cape Split, in the Minas Basin of the Bay of Fundy
cs <- list(longitude=-64.49657, latitude=45.33462)
xy <- lonlat2map(cs, projection="+proj=merc")
map2lonlat(xy)




cleanEx()
nameEx("mapArrows")
### * mapArrows

flush(stderr()); flush(stdout())

### Name: mapArrows
### Title: Add Arrows to a Map
### Aliases: mapArrows

### ** Examples





cleanEx()
nameEx("mapAxis")
### * mapAxis

flush(stderr()); flush(stdout())

### Name: mapAxis
### Title: Add Axis Labels to an Existing Map
### Aliases: mapAxis

### ** Examples





cleanEx()
nameEx("mapContour")
### * mapContour

flush(stderr()); flush(stdout())

### Name: mapContour
### Title: Add Contours on a Existing map
### Aliases: mapContour

### ** Examples

## Not run: 
##D library(oce)
##D data(coastlineWorld)
##D if (requireNamespace("ocedata", quietly=TRUE)) {
##D     data(levitus, package="ocedata")
##D     par(mar=rep(1, 4))
##D     mapPlot(coastlineWorld, projection="+proj=robin", col="lightgray")
##D     mapContour(levitus[['longitude']], levitus[['latitude']], levitus[['SST']])
##D }
## End(Not run)




cleanEx()
nameEx("mapCoordinateSystem")
### * mapCoordinateSystem

flush(stderr()); flush(stdout())

### Name: mapCoordinateSystem
### Title: Draw a coordinate system
### Aliases: mapCoordinateSystem

### ** Examples

## Not run: 
##D library(oce)
##D if (requireNamespace("ocedata", quietly=TRUE)) {
##D     data(coastlineWorldFine, package='ocedata')
##D     HfxLon <- -63.5752
##D     HfxLat <- 44.6488
##D     mapPlot(coastlineWorldFine, proj='+proj=merc',
##D             longitudelim=HfxLon+c(-2,2), latitudelim=HfxLat+c(-2,2),
##D             col='lightgrey')
##D     mapCoordinateSystem(HfxLon, HfxLat, phi=45, length=0.05)
##D    }
## End(Not run)




cleanEx()
nameEx("mapDirectionField")
### * mapDirectionField

flush(stderr()); flush(stdout())

### Name: mapDirectionField
### Title: Add a Direction Field to an Existing Map
### Aliases: mapDirectionField

### ** Examples

## Not run: 
##D library(oce)
##D data(coastlineWorld)
##D par(mar=rep(2, 4))
##D mapPlot(coastlineWorld, longitudelim=c(-120,-55), latitudelim=c(35, 50),
##D         projection="+proj=laea +lat0=40 +lat1=60 +lon_0=-110")
##D lon <- seq(-120, -60, 15)
##D lat <- 45 + seq(-15, 15, 5)
##D lonm <- matrix(expand.grid(lon, lat)[, 1], nrow=length(lon))
##D latm <- matrix(expand.grid(lon, lat)[, 2], nrow=length(lon))
##D ## vectors pointed 45 degrees clockwise from north
##D u <- matrix(1/sqrt(2), nrow=length(lon), ncol=length(lat))
##D v <- matrix(1/sqrt(2), nrow=length(lon), ncol=length(lat))
##D mapDirectionField(lon, lat, u, v, scale=3)
##D mapDirectionField(lonm, latm, 0, 1, scale=3, col='red')
##D # Color code by longitude, using thick lines
##D col <- colormap(lonm)$zcol
##D mapDirectionField(lonm, latm, 1, 0, scale=3, col=col, lwd=2)
## End(Not run)




cleanEx()
nameEx("mapGrid")
### * mapGrid

flush(stderr()); flush(stdout())

### Name: mapGrid
### Title: Add a Longitude and Latitude Grid to an Existing Map
### Aliases: mapGrid

### ** Examples





cleanEx()
nameEx("mapImage")
### * mapImage

flush(stderr()); flush(stdout())

### Name: mapImage
### Title: Add an Image to a Map
### Aliases: mapImage

### ** Examples

## Not run: 
##D library(oce)
##D data(coastlineWorld)
##D data(topoWorld)
##D 
##D ## Northern polar region, with color-coded bathymetry
##D par(mfrow=c(1,1), mar=c(2,2,1,1))
##D cm <- colormap(zlim=c(-5000, 0), col=oceColorsGebco)
##D drawPalette(colormap=cm)
##D mapPlot(coastlineWorld, projection="+proj=stere +lat_0=90",
##D         longitudelim=c(-180,180), latitudelim=c(70,110))
##D mapImage(topoWorld, colormap=cm)
##D mapGrid(15, 15, polarCircle=1, col=gray(0.2))
##D mapPolygon(coastlineWorld[['longitude']], coastlineWorld[['latitude']], col="tan")
## End(Not run)




cleanEx()
nameEx("mapLines")
### * mapLines

flush(stderr()); flush(stdout())

### Name: mapLines
### Title: Add Lines to a Map
### Aliases: mapLines

### ** Examples





cleanEx()
nameEx("mapLongitudeLatitudeXY")
### * mapLongitudeLatitudeXY

flush(stderr()); flush(stdout())

### Name: mapLongitudeLatitudeXY
### Title: Convert From Longitude and Latitude to X and Y
### Aliases: mapLongitudeLatitudeXY

### ** Examples





cleanEx()
nameEx("mapPlot")
### * mapPlot

flush(stderr()); flush(stdout())

### Name: mapPlot
### Title: Draw a Map
### Aliases: mapPlot

### ** Examples


## Not run: 
##D if (utils::packageVersion("sf") != "0.9.8") {
##D     # sf version 0.9-8 has a problem with this projection
##D     png("globe-%03d.png")
##D     lons <- seq(360, 0, -15)
##D     par(mar=rep(0, 4))
##D     for (i in seq_along(lons)) {
##D         p <- paste("+proj=ortho +lat_0=30 +lon_0=", lons[i], sep="")
##D         if (i == 1) {
##D             mapPlot(coastlineCut(coastlineWorld, lons[i]), projection=p, col="gray")
##D             xlim <- par("usr")[1:2]
##D             ylim <- par("usr")[3:4]
##D         } else {
##D             mapPlot(coastlineCut(coastlineWorld, lons[i]), projection=p, col="gray",
##D                     xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
##D         }
##D     }
##D     dev.off()
##D }
## End(Not run)




cleanEx()
nameEx("mapPoints")
### * mapPoints

flush(stderr()); flush(stdout())

### Name: mapPoints
### Title: Add Points to a Map
### Aliases: mapPoints

### ** Examples





cleanEx()
nameEx("mapPolygon")
### * mapPolygon

flush(stderr()); flush(stdout())

### Name: mapPolygon
### Title: Add a Polygon to a Map
### Aliases: mapPolygon

### ** Examples





cleanEx()
nameEx("mapScalebar")
### * mapScalebar

flush(stderr()); flush(stdout())

### Name: mapScalebar
### Title: Add a Scalebar to a Map
### Aliases: mapScalebar

### ** Examples





cleanEx()
nameEx("mapText")
### * mapText

flush(stderr()); flush(stdout())

### Name: mapText
### Title: Add Text to a Map
### Aliases: mapText

### ** Examples





cleanEx()
nameEx("mapTissot")
### * mapTissot

flush(stderr()); flush(stdout())

### Name: mapTissot
### Title: Add Tissot Indicatrices to a Map
### Aliases: mapTissot

### ** Examples





cleanEx()
nameEx("matchBytes")
### * matchBytes

flush(stderr()); flush(stdout())

### Name: matchBytes
### Title: Locate byte sequences in a raw vector
### Aliases: matchBytes

### ** Examples

buf <- as.raw(c(0xa5, 0x11, 0xaa, 0xa5, 0x11, 0x00))
match <- matchBytes(buf, 0xa5, 0x11)
print(buf)
print(match)



cleanEx()
nameEx("matrixSmooth")
### * matrixSmooth

flush(stderr()); flush(stdout())

### Name: matrixSmooth
### Title: Smooth a Matrix
### Aliases: matrixSmooth

### ** Examples

library(oce)
opar <- par(no.readonly = TRUE)
m <- matrix(rep(seq(0, 1, length.out=5), 5), nrow=5, byrow=TRUE)
m[3, 3] <- 2
m1 <- matrixSmooth(m)
m2 <- matrixSmooth(m1)
m3 <- matrixSmooth(m2)
par(mfrow=c(2, 2))
image(m,  col=rainbow(100), zlim=c(0, 4), main="original image")
image(m1, col=rainbow(100), zlim=c(0, 4), main="smoothed 1 time")
image(m2, col=rainbow(100), zlim=c(0, 4), main="smoothed 2 times")
image(m3, col=rainbow(100), zlim=c(0, 4), main="smoothed 3 times")
par(opar)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("moonAngle")
### * moonAngle

flush(stderr()); flush(stdout())

### Name: moonAngle
### Title: Lunar Angle as Function of Space and Time
### Aliases: moonAngle

### ** Examples


library(oce)
par(mfrow=c(3,2))
y <- 2012
m <- 4
days <- 1:3
## Halifax sunrise/sunset (see e.g. https://www.timeanddate.com/worldclock)
rises <- ISOdatetime(y, m, days,c(13,15,16), c(55, 04, 16),0,tz="UTC") + 3 * 3600 # ADT
sets <- ISOdatetime(y, m,days,c(3,4,4), c(42, 15, 45),0,tz="UTC") + 3 * 3600
azrises <- c(69, 75, 82)
azsets <- c(293, 288, 281)
latitude <- 44.65
longitude <- -63.6
for (i in 1:3) {
    t <- ISOdatetime(y, m, days[i],0,0,0,tz="UTC") + seq(0, 24*3600, 3600/4)
    ma <- moonAngle(t, longitude, latitude)
    oce.plot.ts(t, ma$altitude, type='l', mar=c(2, 3, 1, 1), cex=1/2, ylab="Altitude")
    abline(h=0)
    points(rises[i], 0, col='red', pch=3, lwd=2, cex=1.5)
    points(sets[i], 0, col='blue', pch=3, lwd=2, cex=1.5)
    oce.plot.ts(t, ma$azimuth, type='l', mar=c(2, 3, 1, 1), cex=1/2, ylab="Azimuth")
    points(rises[i], -180+azrises[i], col='red', pch=3, lwd=2, cex=1.5)
    points(sets[i], -180+azsets[i], col='blue', pch=3, lwd=2, cex=1.5)
}




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("numberAsHMS")
### * numberAsHMS

flush(stderr()); flush(stdout())

### Name: numberAsHMS
### Title: Convert a Numeric Time to Hour, Minute, and Second
### Aliases: numberAsHMS

### ** Examples

t <- c("0900", "1234")
numberAsHMS(t)




cleanEx()
nameEx("numberAsPOSIXct")
### * numberAsPOSIXct

flush(stderr()); flush(stdout())

### Name: numberAsPOSIXct
### Title: Convert a Numeric Time to a POSIXct Time
### Aliases: numberAsPOSIXct

### ** Examples

numberAsPOSIXct(0)                     # unix time 0
numberAsPOSIXct(1, type="matlab")      # matlab time 1
numberAsPOSIXct(cbind(566, 345615), type="gps") # Canada Day, zero hour UTC
numberAsPOSIXct(cbind(2013, 1), type="yearday") # start of 2013

## Epic time, one hour into Canada Day of year 2018. In computing the
## Julian day, note that this starts at noon.
jd <- julianDay(as.POSIXct("2018-07-01 12:00:00", tz="UTC"))
numberAsPOSIXct(cbind(jd, 1e3 * 1 * 3600), type="epic", tz="UTC")




cleanEx()
nameEx("oce-class")
### * oce-class

flush(stderr()); flush(stdout())

### Name: oce-class
### Title: Base Class for oce Objects
### Aliases: oce-class

### ** Examples

str(new("oce"))




cleanEx()
nameEx("oce.as.raw")
### * oce.as.raw

flush(stderr()); flush(stdout())

### Name: oce.as.raw
### Title: Version of as.raw() that clips data
### Aliases: oce.as.raw

### ** Examples

x <- c(-0.1, 0, 1, 255, 255.1)
data.frame(x, oce.as.raw(x))



cleanEx()
nameEx("oce.contour")
### * oce.contour

flush(stderr()); flush(stdout())

### Name: oce.contour
### Title: Oce Variant of contour
### Aliases: oce.contour oceContour

### ** Examples

library(oce)
data(topoWorld)
## coastline now, and in last glacial maximum
lon <- topoWorld[["longitude"]]
lat <- topoWorld[["latitude"]]
z <- topoWorld[["z"]]
oce.contour(lon, lat, z, levels=0, drawlabels=FALSE)
oce.contour(lon, lat, z, levels=-130, drawlabels=FALSE, col='blue', add=TRUE)



cleanEx()
nameEx("oce.grid")
### * oce.grid

flush(stderr()); flush(stdout())

### Name: oce.grid
### Title: Add a Grid to an Existing Oce Plot
### Aliases: oce.grid

### ** Examples

library(oce)
i <- imagep(volcano)
oce.grid(i, lwd=2)

data(sealevel)
i <- oce.plot.ts(sealevel[["time"]], sealevel[["elevation"]])
oce.grid(i, col='red')

data(ctd)
i <- plotTS(ctd)
oce.grid(i, col='red')

data(adp)
i <- plot(adp, which=1)
oce.grid(i, col='gray', lty=1)

data(echosounder)
i <- plot(echosounder)
oce.grid(i, col='pink', lty=1)




cleanEx()
nameEx("oce.plot.ts")
### * oce.plot.ts

flush(stderr()); flush(stdout())

### Name: oce.plot.ts
### Title: Oce Variant of plot.ts
### Aliases: oce.plot.ts

### ** Examples

library(oce)
t0 <- as.POSIXct("2008-01-01", tz="UTC")
t <- seq(t0, length.out=48, by="30 min")
y <- sin(as.numeric(t - t0) * 2 * pi / (12 * 3600))
oce.plot.ts(t, y, type='l', xaxs='i')
# Show how col, pch and cex get recycled
oce.plot.ts(t, y, type='p', xaxs='i',
            col=1:3, pch=c(rep(1, 6), rep(20, 6)), cex=sqrt(1:6))
# Trimming x; note the narrowing of the y view
oce.plot.ts(t, y, type='p', xlim=c(t[6], t[12]))
# Flip the y axis
oce.plot.ts(t, y, flipy=TRUE)



cleanEx()
nameEx("oceApprox")
### * oceApprox

flush(stderr()); flush(stdout())

### Name: oceApprox
### Title: Interpolate 1D Data with UNESCO or Reiniger-Ross Algorithm
### Aliases: oceApprox oce.approx

### ** Examples

library(oce)
if (require(ocedata)) {
    data(RRprofile)
    zz <- seq(0, 2000, 2)
    plot(RRprofile$temperature, RRprofile$depth, ylim=c(500, 0), xlim=c(2, 11))
    ## Contrast two methods
    a1 <- oce.approx(RRprofile$depth, RRprofile$temperature, zz, "rr")
    a2 <- oce.approx(RRprofile$depth, RRprofile$temperature, zz, "unesco")
    lines(a1, zz)
    lines(a2, zz, col='red')
    legend("bottomright",lwd=1,col=1:2, legend=c("rr","unesco"),cex=3/4)
}



cleanEx()
nameEx("oceAxis")
### * oceAxis

flush(stderr()); flush(stdout())

### Name: oceAxis
### Title: Draw an axis, possibly with decade-style logarithmic scaling
### Aliases: oceAxis

### ** Examples

library(oce)
Ra <- 10^seq(4, 10, 0.1)
Nu <- 0.085 * Ra^(1/3)
plot(Ra, Nu, log="xy", axes=FALSE)
box()
oceAxis(1, logStyle="decade")
oceAxis(2, logStyle="decade")




cleanEx()
nameEx("oceCRS")
### * oceCRS

flush(stderr()); flush(stdout())

### Name: oceCRS
### Title: Coordinate Reference System strings for some oceans
### Aliases: oceCRS

### ** Examples




cleanEx()
nameEx("oceColors9B")
### * oceColors9B

flush(stderr()); flush(stdout())

### Name: oceColors9B
### Title: Create colors in a red-yellow-blue color scheme
### Aliases: oceColors9B oce.colors9B

### ** Examples

library(oce)
imagep(volcano, col=oceColors9B(128),
       zlab="oceColors9B")



cleanEx()
nameEx("oceColorsCDOM")
### * oceColorsCDOM

flush(stderr()); flush(stdout())

### Name: oceColorsCDOM
### Title: Create colors suitable for CDOM fields
### Aliases: oceColorsCDOM oce.colorsCDOM

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsCDOM(128),
       zlab="oceColorsCDOM")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsChlorophyll")
### * oceColorsChlorophyll

flush(stderr()); flush(stdout())

### Name: oceColorsChlorophyll
### Title: Create colors suitable for chlorophyll fields
### Aliases: oceColorsChlorophyll oce.colorsChlorophyll

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsChlorophyll(128),
       zlab="oceColorsChlorophyll")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsClosure")
### * oceColorsClosure

flush(stderr()); flush(stdout())

### Name: oceColorsClosure
### Title: Create color functions
### Aliases: oceColorsClosure

### ** Examples

## Not run: 
##D ## Update oxygen color scheme to latest matplotlib value.
##D library(oce)
##D oxy <- "https://raw.githubusercontent.com/matplotlib/cmocean/master/cmocean/rgb/oxy-rgb.txt"
##D oxyrgb <- read.table(oxy, header=FALSE)
##D oceColorsOxygenUpdated <- oceColorsClosure(oxyrgb)
##D par(mfrow=c(1, 2))
##D m <- matrix(1:256)
##D imagep(m, col=oceColorsOxygen, zlab="oxygen")
##D imagep(m, col=oceColorsOxygenUpdated, zlab="oxygenUpdated")
## End(Not run)



cleanEx()
nameEx("oceColorsDensity")
### * oceColorsDensity

flush(stderr()); flush(stdout())

### Name: oceColorsDensity
### Title: Create colors suitable for density fields
### Aliases: oceColorsDensity oce.colorsDensity

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsDensity(128),
       zlab="oceColorsDensity")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsFreesurface")
### * oceColorsFreesurface

flush(stderr()); flush(stdout())

### Name: oceColorsFreesurface
### Title: Create colors suitable for freesurface fields
### Aliases: oceColorsFreesurface oce.colorsFreesurface

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsFreesurface(128),
       zlab="oceColorsFreesurface")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsGebco")
### * oceColorsGebco

flush(stderr()); flush(stdout())

### Name: oceColorsGebco
### Title: Create colors in a Gebco-like scheme
### Aliases: oceColorsGebco oce.colorsGebco

### ** Examples

library(oce)
imagep(volcano, col=oceColorsGebco(128, region="both"))




cleanEx()
nameEx("oceColorsJet")
### * oceColorsJet

flush(stderr()); flush(stdout())

### Name: oceColorsJet
### Title: Create colors similar to the Matlab Jet scheme
### Aliases: oceColorsJet oce.colorsJet oceColors9A oce.colors9A

### ** Examples

library(oce)
imagep(volcano, col=oceColorsJet, zlab="oceColorsJet")



cleanEx()
nameEx("oceColorsOxygen")
### * oceColorsOxygen

flush(stderr()); flush(stdout())

### Name: oceColorsOxygen
### Title: Create colors suitable for oxygen fields
### Aliases: oceColorsOxygen oce.colorsOxygen

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsOxygen(128),
       zlab="oceColorsOxygen")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsPAR")
### * oceColorsPAR

flush(stderr()); flush(stdout())

### Name: oceColorsPAR
### Title: Create colors suitable for PAR fields
### Aliases: oceColorsPAR oce.colorsPAR

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsPAR(128),
       zlab="oceColorsPAR")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsPhase")
### * oceColorsPhase

flush(stderr()); flush(stdout())

### Name: oceColorsPhase
### Title: Create colors suitable for phase fields
### Aliases: oceColorsPhase oce.colorsPhase

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsPhase(128),
       zlab="oceColorsPhase")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsSalinity")
### * oceColorsSalinity

flush(stderr()); flush(stdout())

### Name: oceColorsSalinity
### Title: Create colors suitable for salinity fields
### Aliases: oceColorsSalinity oce.colorsSalinity

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsSalinity(128),
       zlab="oceColorsSalinity")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsTemperature")
### * oceColorsTemperature

flush(stderr()); flush(stdout())

### Name: oceColorsTemperature
### Title: Create colors suitable for temperature fields
### Aliases: oceColorsTemperature oce.colorsTemperature

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsTemperature(128),
       zlab="oceColorsTemperature")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsTurbidity")
### * oceColorsTurbidity

flush(stderr()); flush(stdout())

### Name: oceColorsTurbidity
### Title: Create colors suitable for turbidity fields
### Aliases: oceColorsTurbidity oce.colorsTurbidity

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsTurbidity(128),
       zlab="oceColorsTurbidity")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsTurbo")
### * oceColorsTurbo

flush(stderr()); flush(stdout())

### Name: oceColorsTurbo
### Title: Create colors similar to the google turbo scheme
### Aliases: oceColorsTurbo oce.colorsTurbo

### ** Examples

library(oce)
imagep(volcano, col=oceColorsTurbo(128),
       zlab="oceColorsTurbo")



cleanEx()
nameEx("oceColorsTwo")
### * oceColorsTwo

flush(stderr()); flush(stdout())

### Name: oceColorsTwo
### Title: Create two-color palette
### Aliases: oceColorsTwo oce.colorsTwo

### ** Examples

library(oce)
imagep(volcano-mean(range(volcano)), col=oceColorsTwo(128),
       zlim="symmetric", zlab="oceColorsTwo")



cleanEx()
nameEx("oceColorsVelocity")
### * oceColorsVelocity

flush(stderr()); flush(stdout())

### Name: oceColorsVelocity
### Title: Create colors suitable for velocity fields
### Aliases: oceColorsVelocity oce.colorsVelocity

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsVelocity(128),
       zlab="oceColorsVelocity")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceColorsViridis")
### * oceColorsViridis

flush(stderr()); flush(stdout())

### Name: oceColorsViridis
### Title: Create colors similar to the matlab Viridis scheme
### Aliases: oceColorsViridis oce.colorsViridis

### ** Examples

library(oce)
# Example 1: oceColorsViridis
imagep(volcano, col=oceColorsViridis(128),
       zlab="oceColorsViridis")
# Example 2: viridis::viridis
## Not run: 
##D imagep(volcano, col=viridis::viridis,
##D        zlab="viridis::viridis")
## End(Not run)




cleanEx()
nameEx("oceColorsVorticity")
### * oceColorsVorticity

flush(stderr()); flush(stdout())

### Name: oceColorsVorticity
### Title: Create colors suitable for vorticity fields
### Aliases: oceColorsVorticity oce.colorsVorticity

### ** Examples

library(oce)

# Example 1
imagep(volcano, col=oceColorsVorticity(128),
       zlab="oceColorsVorticity")
## Not run: 
##D # Example 2 (requires the cmocean package)
##D imagep(volcano, col=cmocean::cmocean("matter"),
##D        zlab='cmocean::cmocean("matter")')
## End(Not run)

## Not run: 
##D # Example 3 (requires the viridis package)
##D imagep(volcano, col=viridis::inferno,
##D        zlab='viridis::inferno')
## End(Not run)




cleanEx()
nameEx("oceConvolve")
### * oceConvolve

flush(stderr()); flush(stdout())

### Name: oceConvolve
### Title: Convolve two time series
### Aliases: oceConvolve oce.convolve

### ** Examples

library(oce)
t <- 0:1027
n <- length(t)
signal <- ifelse(sin(t * 2 * pi / 128) > 0, 1, 0)
tau <- 10
filter <- exp(-seq(5*tau, 0) / tau)
filter <- filter / sum(filter)
observation <- oce.convolve(signal, filter)
plot(t, signal, type='l')
lines(t, observation, lty='dotted')



cleanEx()
nameEx("oceDebug")
### * oceDebug

flush(stderr()); flush(stdout())

### Name: oceDebug
### Title: Print a debugging message
### Aliases: oceDebug oce.debug

### ** Examples

oceDebug(debug=1, "Example", 1, "Plain text")
oceDebug(debug=1, "Example", 2, "Bold", style="bold")
oceDebug(debug=1, "Example", 3, "Italic", style="italic")
oceDebug(debug=1, "Example", 4, "Red", style="red")
oceDebug(debug=1, "Example", 5, "Green", style="green")
oceDebug(debug=1, "Example", 6, "Blue", style="blue")
mycyan <- function(...) paste("\033[36m", paste(..., sep=" "), "\033[0m", sep="")
oceDebug(debug=1, "Example", 7, "User-set cyan", style=mycyan)



cleanEx()
nameEx("oceEdit")
### * oceEdit

flush(stderr()); flush(stdout())

### Name: oceEdit
### Title: Edit an Oce Object
### Aliases: oceEdit oce.edit

### ** Examples

library(oce)
data(ctd)
ctd2 <- oceEdit(ctd, item="latitude", value=47.8879,
               reason="illustration", person="Dan Kelley")
ctd3 <- oceEdit(ctd,action="x@data$pressure<-x@data$pressure-1")



cleanEx()
nameEx("oceFilter")
### * oceFilter

flush(stderr()); flush(stdout())

### Name: oceFilter
### Title: Filter a Time Series
### Aliases: oceFilter oce.filter

### ** Examples

library(oce)
par(mar=c(4, 4, 1, 1))
b <- rep(1, 5)/5
a <- 1
x <- seq(0, 10)
y <- ifelse(x == 5, 1, 0)
f1 <- oceFilter(y, a, b)
plot(x, y, ylim=c(-0, 1.5), pch="o", type='b')
points(x, f1, pch="x", col="red")

# remove the phase lag
f2 <- oceFilter(y, a, b, TRUE)
points(x, f2, pch="+", col="blue")

legend("topleft", col=c("black","red","blue"), pch=c("o","x","+"),
       legend=c("data","normal filter", "zero-phase filter"))
mtext("note that normal filter rolls off at end")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ocePmatch")
### * ocePmatch

flush(stderr()); flush(stdout())

### Name: ocePmatch
### Title: Partial matching of strings or numbers
### Aliases: ocePmatch oce.pmatch

### ** Examples

library(oce)
oce.pmatch(c("s", "at", "te"), list(salinity=1, temperature=3.1))



cleanEx()
nameEx("oceRenameData")
### * oceRenameData

flush(stderr()); flush(stdout())

### Name: oceRenameData
### Title: Rename Something in an oce data Slot
### Aliases: oceRenameData

### ** Examples

library(oce)
data(ctd)
CTD <- oceRenameData(ctd, "salinity", "SALT")
expect_equal(ctd[["salinity"]], CTD[["SALT"]])
expect_equal(ctd[["sal00"]], CTD[["SALT"]])




cleanEx()
nameEx("oceSetData")
### * oceSetData

flush(stderr()); flush(stdout())

### Name: oceSetData
### Title: Set Something in an oce data Slot
### Aliases: oceSetData

### ** Examples

data(ctd)
Tf <- swTFreeze(ctd)
ctd <- oceSetData(ctd, "freezing", Tf,
    unit=list(unit=expression(degree*C), scale="ITS-90"))
plotProfile(ctd, "freezing")




cleanEx()
nameEx("oceSetMetadata")
### * oceSetMetadata

flush(stderr()); flush(stdout())

### Name: oceSetMetadata
### Title: Set Something in an oce metadata Slot
### Aliases: oceSetMetadata

### ** Examples

# Add an estimate of MLD (mixed layer depth) to a ctd object
library(oce)
data(ctd)
ctdWithMLD <- oceSetMetadata(ctd, "MLD", 3)
ctdWithMLD[["MLD"]] # 3




cleanEx()
nameEx("oceSmooth")
### * oceSmooth

flush(stderr()); flush(stdout())

### Name: oceSmooth
### Title: Smooth an Oce Object
### Aliases: oceSmooth oce.smooth

### ** Examples

library(oce)
data(ctd)
d <- oce.smooth(ctd)
plot(d)



cleanEx()
nameEx("oceSpectrum")
### * oceSpectrum

flush(stderr()); flush(stdout())

### Name: oceSpectrum
### Title: Wrapper to give normalized spectrum
### Aliases: oceSpectrum oce.spectrum

### ** Examples

  x <- rnorm(1e3)
  s <- spectrum(x, plot=FALSE)
  ss <- oce.spectrum(x, plot=FALSE)
  cat("variance of x=", var(x), "\n")
  cat("integral of     spectrum=", sum(s$spec)*diff(s$freq[1:2]), "\n")
  cat("integral of oce.spectrum=", sum(ss$spec)*diff(ss$freq[1:2]), "\n")



cleanEx()
nameEx("plot-adp-method")
### * plot-adp-method

flush(stderr()); flush(stdout())

### Name: plot,adp-method
### Title: Plot an adp Object
### Aliases: plot,adp-method plot.adp

### ** Examples

library(oce)
data(adp)
plot(adp, which=1:3)
plot(adp, which='temperature', tformat='%H:%M')




cleanEx()
nameEx("plot-adv-method")
### * plot-adv-method

flush(stderr()); flush(stdout())

### Name: plot,adv-method
### Title: Plot an adv Object
### Aliases: plot,adv-method plot.adv

### ** Examples

library(oce)
data(adv)
plot(adv)




cleanEx()
nameEx("plot-amsr-method")
### * plot-amsr-method

flush(stderr()); flush(stdout())

### Name: plot,amsr-method
### Title: Plot an amsr Object
### Aliases: plot,amsr-method plot.amsr both xlim and ylim to compute
###   zrange

### ** Examples

library(oce)
data(coastlineWorld)
data(amsr) # see ?amsr for how to read and composite such objects

# Example 1: plot with default colour scheme, oceColorsTemperature()
plot(amsr, "SST")
lines(coastlineWorld[['longitude']], coastlineWorld[['latitude']])

# Example 2: 'turbo' colour scheme
plot(amsr, "SST", col=oceColorsTurbo)
lines(coastlineWorld[['longitude']], coastlineWorld[['latitude']])




cleanEx()
nameEx("plot-argo-method")
### * plot-argo-method

flush(stderr()); flush(stdout())

### Name: plot,argo-method
### Title: Plot an argo Object
### Aliases: plot,argo-method plot.argo

### ** Examples

library(oce)
data(argo)
tc <- cut(argo[["time"]], "year")
# Example 1: plot map, which reveals float trajectory.
plot(argo, pch=as.integer(tc))
year <- substr(levels(tc), 1, 4)
data(topoWorld)
contour(topoWorld[['longitude']], topoWorld[['latitude']],
        topoWorld[['z']], add=TRUE)
legend("bottomleft", pch=seq_along(year), legend=year, bg="white", cex=3/4)

# Example 2: plot map, TS, T(z) and S(z). Note the use
# of handleFlags(), to skip over questionable data.
plot(handleFlags(argo), which=c(1, 4, 6, 5))




cleanEx()
nameEx("plot-cm-method")
### * plot-cm-method

flush(stderr()); flush(stdout())

### Name: plot,cm-method
### Title: Plot a cm Object
### Aliases: plot,cm-method plot.cm

### ** Examples

  library(oce)
  data(cm)
  summary(cm)
  plot(cm)




cleanEx()
nameEx("plot-coastline-method")
### * plot-coastline-method

flush(stderr()); flush(stdout())

### Name: plot,coastline-method
### Title: Plot a coastline Object
### Aliases: plot,coastline-method plot.coastline

### ** Examples





cleanEx()
nameEx("plot-ctd-method")
### * plot-ctd-method

flush(stderr()); flush(stdout())

### Name: plot,ctd-method
### Title: Plot a ctd Object
### Aliases: plot,ctd-method plot.ctd

### ** Examples

## 1. simple plot
library(oce)
data(ctd)
plot(ctd)

## 2. how to customize depth contours
par(mfrow=c(1,2))
data(section)
stn <- section[["station", 105]]
plot(stn, which='map', drawIsobaths=TRUE)
plot(stn, which='map')
data(topoWorld)
tlon <- topoWorld[["longitude"]]
tlat <- topoWorld[["latitude"]]
tdep <- -topoWorld[["z"]]
contour(tlon, tlat, tdep, drawlabels=FALSE,
        levels=seq(1000,6000,1000), col='lightblue', add=TRUE)
contour(tlon, tlat, tdep, vfont=c("sans serif", "bold"),
        levels=stn[['waterDepth']], col='red', lwd=2, add=TRUE)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plot-echosounder-method")
### * plot-echosounder-method

flush(stderr()); flush(stdout())

### Name: plot,echosounder-method
### Title: Plot an echosounder Object
### Aliases: plot,echosounder-method plot.echosounder

### ** Examples

library(oce)
data(echosounder)
plot(echosounder, drawBottom=TRUE)




cleanEx()
nameEx("plot-lisst-method")
### * plot-lisst-method

flush(stderr()); flush(stdout())

### Name: plot,lisst-method
### Title: Plot a lisst Object
### Aliases: plot,lisst-method plot.lisst

### ** Examples

library(oce)
data(lisst)
plot(lisst)




cleanEx()
nameEx("plot-met-method")
### * plot-met-method

flush(stderr()); flush(stdout())

### Name: plot,met-method
### Title: Plot a met Object
### Aliases: plot,met-method plot.met

### ** Examples

library(oce)
data(met)
plot(met, which=3:4)

## Wind speed and direction during Hurricane Juan
## Compare with the final figure in a white paper by Chris Fogarty
## (available at http://www.novaweather.net/Hurricane_Juan_files/McNabs_plot.pdf
## downloaded 2017-01-02).
library(oce)
data(met)
t0 <- as.POSIXct("2003-09-29 04:00:00", tz="UTC")
dt <- 12 * 3600
juan <- subset(met, t0 - dt <= time & time <= t0 + dt)
par(mfrow=c(2,1))
plot(juan, which=5)
abline(v=t0)
plot(juan, which=6)
abline(v=t0)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plot-oce-method")
### * plot-oce-method

flush(stderr()); flush(stdout())

### Name: plot,oce-method
### Title: Plot an oce Object
### Aliases: plot,oce-method plot.oce

### ** Examples

library(oce)
o <- new("oce")
o <- oceSetData(o, 'x', rnorm(10))
o <- oceSetData(o, 'y', rnorm(10))
o <- oceSetData(o, 'z', rnorm(10))
plot(o)



cleanEx()
nameEx("plot-rsk-method")
### * plot-rsk-method

flush(stderr()); flush(stdout())

### Name: plot,rsk-method
### Title: Plot a rsk Object
### Aliases: plot,rsk-method plot.rsk

### ** Examples

library(oce)
data(rsk)
plot(rsk) # default timeseries plot of all data fields

## A multipanel plot of just pressure and temperature with ylim
par(mfrow=c(2, 1))
plot(rsk, which="pressure", ylim=c(10, 30))
plot(rsk, which="temperature", ylim=c(2, 4))




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plot-sealevel-method")
### * plot-sealevel-method

flush(stderr()); flush(stdout())

### Name: plot,sealevel-method
### Title: Plot a sealevel Object
### Aliases: plot,sealevel-method plot.sealevel

### ** Examples

library(oce)
data(sealevel)
## local Halifax time is UTC + 4h
juan <- as.POSIXct("2003-09-29 00:15:00", tz="UTC")+4*3600
plot(sealevel, which=1, xlim=juan+86400*c(-7, 7))
abline(v=juan, col='red')




cleanEx()
nameEx("plot-section-method")
### * plot-section-method

flush(stderr()); flush(stdout())

### Name: plot,section-method
### Title: Plot a section Object
### Aliases: plot,section-method plot.section

### ** Examples

library(oce)
data(section)
sg <- sectionGrid(section)

## 1. start of section, default fields.
plot(head(section))

## 2. Gulf Stream
GS <- subset(section, 109<=stationId&stationId<=129)
GSg <- sectionGrid(GS, p=seq(0, 2000, 100))
plot(GSg, which=c(1, 99), map.ylim=c(34, 42))
par(mfrow=c(2, 1))
plot(GS, which=1, ylim=c(2000, 0), ztype='points',
     zbreaks=seq(0,30,2), pch=20, cex=3)
plot(GSg, which=1, ztype='image', zbreaks=seq(0,30,2))

par(mfrow=c(1, 1))

## 3. Image, with colored dots to indicate grid-data mismatch.
## Not run: 
##D plot(GSg, which=1, ztype='image')
##D T <- GS[['temperature']]
##D col <- oceColorsViridis(100)[rescale(T, rlow=1, rhigh=100)]
##D points(GS[['distance']],GS[['depth']],pch=20,cex=3,col='white')
##D points(GS[['distance']],GS[['depth']],pch=20,cex=2.5,col=col)
## End(Not run)


## Not run: 
##D ## 4. Image of Absolute Salinity, with 4-minute bathymetry
##D ## It's easy to calculate the desired area for the bathymetry,
##D ## but for brevity we'll hard-code it. Note that download.topo()
##D ## caches the file locally.
##D f <- download.topo(west=-80, east=0, south=35, north=40, resolution=4)
##D t <- read.topo(f)
##D plot(section, which="SA", xtype="longitude", ztype="image", showBottom=t)
## End(Not run)

## Not run: 
##D ## 5. Temperature with salinity added in red
##D s <- plot(section, which="temperature")
##D distance <- s[["distance", "byStation"]]
##D depth <- s[["station", 1]][["depth"]]
##D salinity <- matrix(s[["salinity"]], byrow=TRUE, nrow=length(s[["station"]]))
##D contour(distance, depth, salinity, col=2, add=TRUE)
## End(Not run)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plot-tidem-method")
### * plot-tidem-method

flush(stderr()); flush(stdout())

### Name: plot,tidem-method
### Title: Plot a tidem Object
### Aliases: plot,tidem-method plot.tidem

### ** Examples

## Not run: 
##D library(oce)
##D data(sealevel)
##D tide <- tidem(sealevel)
##D plot(tide)
## End(Not run)




cleanEx()
nameEx("plot-topo-method")
### * plot-topo-method

flush(stderr()); flush(stdout())

### Name: plot,topo-method
### Title: Plot a topo Object
### Aliases: plot,topo-method plot.topo

### ** Examples

library(oce)
data(topoWorld)
plot(topoWorld, clongitude=-60, clatitude=45, span=10000)




cleanEx()
nameEx("plot-windrose-method")
### * plot-windrose-method

flush(stderr()); flush(stdout())

### Name: plot,windrose-method
### Title: Plot a windrose Object
### Aliases: plot,windrose-method plot.windrose

### ** Examples

library(oce)
opar <- par(no.readonly = TRUE)
xcomp <- rnorm(360) + 1
ycomp <- rnorm(360)
wr <- as.windrose(xcomp, ycomp)
par(mfrow=c(1, 2))
plot(wr)
plot(wr, "fivenum")
par(opar)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plot-xbt-method")
### * plot-xbt-method

flush(stderr()); flush(stdout())

### Name: plot,xbt-method
### Title: Plot an xbt Object
### Aliases: plot,xbt-method plot.xbt

### ** Examples

library(oce)
data(xbt)
summary(xbt)
plot(xbt)




cleanEx()
nameEx("plotInset")
### * plotInset

flush(stderr()); flush(stdout())

### Name: plotInset
### Title: Plot an Inset Diagram
### Aliases: plotInset

### ** Examples

library(oce)
## power law in linear and log form
x <- 1:10
y <- x^2
plot(x, y, log='xy',type='l')
plotInset(3, 1, 10, 8,
          expr=plot(x,y,type='l',cex.axis=3/4,mgp=c(3/2, 1/2, 0)),
          mar=c(2.5, 2.5, 1, 1))

## CTD data with location
data(ctd)
plot(ctd, which="TS")
plotInset(29.9, 2.7, 31, 10,
          expr=plot(ctd, which='map',
          coastline="coastlineWorld",
          span=5000, mar=NULL, cex.axis=3/4))



cleanEx()
nameEx("plotPolar")
### * plotPolar

flush(stderr()); flush(stdout())

### Name: plotPolar
### Title: Draw a Polar Plot
### Aliases: plotPolar

### ** Examples

library(oce)
r <- rnorm(50, mean=2, sd=0.1)
theta <- runif(50, 0, 360)
plotPolar(r, theta)



cleanEx()
nameEx("plotProfile")
### * plotProfile

flush(stderr()); flush(stdout())

### Name: plotProfile
### Title: Plot a CTD Profile
### Aliases: plotProfile

### ** Examples

library(oce)
data(ctd)
plotProfile(ctd, xtype="temperature")




cleanEx()
nameEx("plotScan")
### * plotScan

flush(stderr()); flush(stdout())

### Name: plotScan
### Title: Plot CTD data in a Low-Level Fashion
### Aliases: plotScan

### ** Examples

library(oce)
data(ctdRaw)
plotScan(ctdRaw)
abline(v=c(130, 350), col='red') # useful for ctdTrim()




cleanEx()
nameEx("plotSticks")
### * plotSticks

flush(stderr()); flush(stdout())

### Name: plotSticks
### Title: Draw a Stick Plot
### Aliases: plotSticks

### ** Examples

library(oce)

# Flow from a point source
n <- 16
x <- rep(0, n)
y <- rep(0, n)
theta <- seq(0, 2*pi, length.out=n)
u <- sin(theta)
v <- cos(theta)
plotSticks(x, y, u, v, xlim=c(-2, 2), ylim=c(-2, 2))
rm(n, x, y, theta, u, v)

# Oceanographic example
data(met)
t <- met[["time"]]
u <- met[["u"]]
v <- met[["v"]]
p <- met[["pressure"]]
oce.plot.ts(t, p)
plotSticks(t, 99, u, v, yscale=25, add=TRUE)



cleanEx()
nameEx("plotTS")
### * plotTS

flush(stderr()); flush(stdout())

### Name: plotTS
### Title: Plot Temperature-Salinity Diagram
### Aliases: plotTS

### ** Examples

library(oce)
data(ctd)
plotTS(ctd)




cleanEx()
nameEx("plotTaylor")
### * plotTaylor

flush(stderr()); flush(stdout())

### Name: plotTaylor
### Title: Plot a Model-data Comparison Diagram
### Aliases: plotTaylor

### ** Examples

library(oce)
data(sealevel)
x <- sealevel[["elevation"]]
M2 <- predict(tidem(sealevel, constituents="M2"))
S2 <- predict(tidem(sealevel, constituents=c("S2")))
plotTaylor(x, cbind(M2, S2))



cleanEx()
nameEx("predict.tidem")
### * predict.tidem

flush(stderr()); flush(stdout())

### Name: predict.tidem
### Title: Predict a Tidal Signal
### Aliases: predict.tidem

### ** Examples


## Not run: 
##D library(oce)
##D # 1. tidal anomaly
##D data(sealevelTuktoyaktuk)
##D time <- sealevelTuktoyaktuk[["time"]]
##D elevation <- sealevelTuktoyaktuk[["elevation"]]
##D oce.plot.ts(time, elevation, type='l', ylab="Height [m]", ylim=c(-2, 6))
##D tide <- tidem(sealevelTuktoyaktuk)
##D lines(time, elevation - predict(tide), col="red")
##D abline(h=0, col="red")
##D 
##D # 2. prediction at specified times
##D data(sealevel)
##D m <- tidem(sealevel)
##D ## Check fit over 2 days (interpolating to finer timescale)
##D look <- 1:48
##D time <- sealevel[["time"]]
##D elevation <- sealevel[["elevation"]]
##D oce.plot.ts(time[look], elevation[look])
##D # 360s = 10 minute timescale
##D t <- seq(from=time[1], to=time[max(look)], by=360)
##D lines(t, predict(m, newdata=t), col='red')
##D legend("topright", col=c("black","red"),
##D legend=c("data","model"),lwd=1)
## End(Not run)




cleanEx()
nameEx("preferAdjusted")
### * preferAdjusted

flush(stderr()); flush(stdout())

### Name: preferAdjusted
### Title: Set Preference for Adjusted Values
### Aliases: preferAdjusted

### ** Examples

library(oce)
data(argo)
argoAdjusted <- preferAdjusted(argo)
all.equal(argo[["salinityAdjusted"]], argoAdjusted[["salinity"]])
all.equal(argo[["salinityFlagsAdjusted"]], argoAdjusted[["salinityFlags"]])
all.equal(argo[["salinityUnitsAdjusted"]], argoAdjusted[["salinityUnits"]])




cleanEx()
nameEx("presentTime")
### * presentTime

flush(stderr()); flush(stdout())

### Name: presentTime
### Title: Get the present time, in a stated timezone
### Aliases: presentTime

### ** Examples

presentTime() # UTC
presentTime("") # the local timezone




cleanEx()
nameEx("prettyPosition")
### * prettyPosition

flush(stderr()); flush(stdout())

### Name: prettyPosition
### Title: Pretty lat/lon in deg, min, sec
### Aliases: prettyPosition

### ** Examples

library(oce)
formatPosition(prettyPosition(10+1:10/60+2.8/3600))



cleanEx()
nameEx("processingLog-set")
### * processingLog-set

flush(stderr()); flush(stdout())

### Name: processingLog<-
### Title: Add an item to a processing log (in place)
### Aliases: processingLog<-

### ** Examples

data(ctd)
processingLogShow(ctd)
processingLog(ctd) <- "test"
processingLogShow(ctd)




cleanEx()
nameEx("pwelch")
### * pwelch

flush(stderr()); flush(stdout())

### Name: pwelch
### Title: Welch periodogram
### Aliases: pwelch

### ** Examples

library(oce)
Fs <- 1000
t <- seq(0, 0.296, 1/Fs)
x <- cos(2 * pi * t * 200) + rnorm(n=length(t))
X <- ts(x, frequency=Fs)
s <- spectrum(X, spans=c(3,2), main="random + 200 Hz", log='no')
w <- pwelch(X, plot=FALSE)
lines(w$freq, w$spec, col="red")
w2 <- pwelch(X, nfft=75, plot=FALSE)
lines(w2$freq, w2$spec, col='green')
abline(v=200, col="blue", lty="dotted")
cat("Checking spectral levels with Parseval's theorem:\n")
cat("var(x)                              = ", var(x), "\n")
cat("2 * sum(s$spec) * diff(s$freq[1:2]) = ", 2 * sum(s$spec) * diff(s$freq[1:2]), "\n")
cat("sum(w$spec) * diff(s$freq[1:2])     = ", sum(w$spec) * diff(w$freq[1:2]), "\n")
cat("sum(w2$spec) * diff(s$freq[1:2])    = ", sum(w2$spec) * diff(w2$freq[1:2]), "\n")
## co2
par(mar=c(3,3,2,1), mgp=c(2,0.7,0))
s <- spectrum(co2, plot=FALSE)
plot(log10(s$freq), s$spec * s$freq,
     xlab=expression(log[10]*Frequency), ylab="Power*Frequency", type='l')
title("Variance-preserving spectrum")
pw <- pwelch(co2, nfft=256, plot=FALSE)
lines(log10(pw$freq), pw$spec * pw$freq, col='red')



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("rangeLimit")
### * rangeLimit

flush(stderr()); flush(stdout())

### Name: rangeLimit
### Title: Substitute NA for data outside a range
### Aliases: rangeLimit

### ** Examples


ten.to.twenty <- rangeLimit(1:100, 10, 20)



cleanEx()
nameEx("read.adp.ad2cp")
### * read.adp.ad2cp

flush(stderr()); flush(stdout())

### Name: read.adp.ad2cp
### Title: Read an AD2CP File
### Aliases: read.adp.ad2cp

### ** Examples

## Not run: 
##D d <- read.adp.ad2cp("~/test.ad2cp", to=100) # or read.oce()
## End(Not run)




cleanEx()
nameEx("read.adp.rdi")
### * read.adp.rdi

flush(stderr()); flush(stdout())

### Name: read.adp.rdi
### Title: Read a Teledyne/RDI ADP File
### Aliases: read.adp.rdi

### ** Examples

adp <- read.adp.rdi(system.file("extdata", "adp_rdi.000", package="oce"))
summary(adp)




cleanEx()
nameEx("read.adv")
### * read.adv

flush(stderr()); flush(stdout())

### Name: read.adv
### Title: Read an ADV data file
### Aliases: read.adv

### ** Examples

## Not run: 
##D library(oce)
##D # A nortek Vector file
##D d <- read.oce("/data/archive/sleiwex/2008/moorings/m05/adv/nortek_1943/raw/adv_nortek_1943.vec",
##D               from=as.POSIXct("2008-06-26 00:00:00", tz="UTC"),
##D               to=as.POSIXct("2008-06-26 00:00:10", tz="UTC"))
##D plot(d, which=c(1:3,15))
## End(Not run)



cleanEx()
nameEx("read.adv.nortek")
### * read.adv.nortek

flush(stderr()); flush(stdout())

### Name: read.adv.nortek
### Title: Read an ADV data file
### Aliases: read.adv.nortek

### ** Examples

## Not run: 
##D library(oce)
##D # A nortek Vector file
##D d <- read.oce("/data/archive/sleiwex/2008/moorings/m05/adv/nortek_1943/raw/adv_nortek_1943.vec",
##D               from=as.POSIXct("2008-06-26 00:00:00", tz="UTC"),
##D               to=as.POSIXct("2008-06-26 00:00:10", tz="UTC"))
##D plot(d, which=c(1:3,15))
## End(Not run)



cleanEx()
nameEx("read.adv.sontek.adr")
### * read.adv.sontek.adr

flush(stderr()); flush(stdout())

### Name: read.adv.sontek.adr
### Title: Read an ADV data file
### Aliases: read.adv.sontek.adr

### ** Examples

## Not run: 
##D library(oce)
##D # A nortek Vector file
##D d <- read.oce("/data/archive/sleiwex/2008/moorings/m05/adv/nortek_1943/raw/adv_nortek_1943.vec",
##D               from=as.POSIXct("2008-06-26 00:00:00", tz="UTC"),
##D               to=as.POSIXct("2008-06-26 00:00:10", tz="UTC"))
##D plot(d, which=c(1:3,15))
## End(Not run)



cleanEx()
nameEx("read.adv.sontek.serial")
### * read.adv.sontek.serial

flush(stderr()); flush(stdout())

### Name: read.adv.sontek.serial
### Title: Read an ADV data file
### Aliases: read.adv.sontek.serial

### ** Examples

## Not run: 
##D library(oce)
##D # A nortek Vector file
##D d <- read.oce("/data/archive/sleiwex/2008/moorings/m05/adv/nortek_1943/raw/adv_nortek_1943.vec",
##D               from=as.POSIXct("2008-06-26 00:00:00", tz="UTC"),
##D               to=as.POSIXct("2008-06-26 00:00:10", tz="UTC"))
##D plot(d, which=c(1:3,15))
## End(Not run)



cleanEx()
nameEx("read.adv.sontek.text")
### * read.adv.sontek.text

flush(stderr()); flush(stdout())

### Name: read.adv.sontek.text
### Title: Read an ADV data file
### Aliases: read.adv.sontek.text

### ** Examples

## Not run: 
##D library(oce)
##D # A nortek Vector file
##D d <- read.oce("/data/archive/sleiwex/2008/moorings/m05/adv/nortek_1943/raw/adv_nortek_1943.vec",
##D               from=as.POSIXct("2008-06-26 00:00:00", tz="UTC"),
##D               to=as.POSIXct("2008-06-26 00:00:10", tz="UTC"))
##D plot(d, which=c(1:3,15))
## End(Not run)



cleanEx()
nameEx("read.argo")
### * read.argo

flush(stderr()); flush(stdout())

### Name: read.argo
### Title: Read an Argo Data File
### Aliases: read.argo

### ** Examples

## Not run: 
##D ## Example 1: read from a local file
##D library(oce)
##D d <- read.argo("/data/OAR/6900388_prof.nc")
##D summary(d)
##D plot(d)
##D 
##D ## Example 2: construct URL for download (brittle)
##D id <- "6900388"
##D url <- "https://www.usgodae.org/ftp/outgoing/argo"
##D if (!length(list.files(pattern="argo_index.txt")))
##D     download.file(paste(url, "ar_index_global_meta.txt", sep="/"), "argo_index.txt")
##D index <- readLines("argo_index.txt")
##D line <- grep(id, index)
##D if (0 == length(line)) stop("id ", id, " not found")
##D if (1 < length(line)) stop("id ", id, " found multiple times")
##D dac <- strsplit(index[line], "/")[[1]][1]
##D profile <- paste(id, "_prof.nc", sep="")
##D float <- paste(url, "dac", dac, id, profile, sep="/")
##D download.file(float, profile)
##D argo <- read.argo(profile)
##D summary(argo)
## End(Not run)





cleanEx()
nameEx("read.cm")
### * read.cm

flush(stderr()); flush(stdout())

### Name: read.cm
### Title: Read a CM file
### Aliases: read.cm

### ** Examples

## Not run: 
##D   library(oce)
##D   cm <- read.oce("cm_interocean_0811786.s4a.tab")
##D   summary(cm)
##D   plot(cm)
## End(Not run)




cleanEx()
nameEx("read.ctd.sbe")
### * read.ctd.sbe

flush(stderr()); flush(stdout())

### Name: read.ctd.sbe
### Title: Read a Seabird CTD File
### Aliases: read.ctd.sbe

### ** Examples

f <- system.file("extdata", "ctd.cnv", package="oce")
## Read the file in the normal way
d <- read.ctd(f)
## Read an imaginary file, in which salinity is named 'salt'
d <- read.ctd(f, columns=list(
  salinity=list(name="salt", unit=list(unit=expression(), scale="PSS-78"))))




cleanEx()
nameEx("read.g1sst")
### * read.g1sst

flush(stderr()); flush(stdout())

### Name: read.g1sst
### Title: Read a G1SST file
### Aliases: read.g1sst

### ** Examples

## Not run: 
##D # Construct query, making it easier to understand and modify.
##D day <- "2016-01-02"
##D lon0 <- -66.5
##D lon1 <- -64.0
##D lat0 <- 44
##D lat1 <- 46
##D source <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/",
##D                 "jplG1SST.nc?",
##D                 "SST%5B(", day, "T12:00:00Z)",
##D                 "%5D%5B(", lat0, "):(", lat1, ")",
##D                 "%5D%5B(", lon0, "):(", lon1, ")",
##D                 "%5D", sep="")
##D if (!length(list.files(pattern="^a.nc$")))
##D     download.file(source, "a.nc")
##D d <- read.g1sst("a.nc")
##D plot(d, "SST", col=oceColorsTemperature)
##D if (requireNamespace("ocedata", quietly=TRUE)) {
##D     data(coastlineWorldFine, package="ocedata")
##D     lines(coastlineWorldFine[['longitude']],coastlineWorldFine[['latitude']])
##D }
## End(Not run)




cleanEx()
nameEx("read.index")
### * read.index

flush(stderr()); flush(stdout())

### Name: read.index
### Title: Read a NOAA ocean index file
### Aliases: read.index

### ** Examples


## Not run: 
##D library(oce)
##D par(mfrow=c(2, 1))
##D # 1. AO, Arctic oscillation
##D #  Note that data used to be at https://www.esrl.noaa.gov/psd/data/correlation/ao.data
##D ao <- read.index("https://psl.noaa.gov/data/correlation/ao.data")
##D aorecent <- subset(ao, t > as.POSIXct("2000-01-01"))
##D oce.plot.ts(aorecent$t, aorecent$index)
##D # 2. SOI, probably more up-to-date then data(soi, package="ocedata")
##D soi <- read.index("https://www.cgd.ucar.edu/cas/catalog/climind/SOI.signal.ascii")
##D soirecent <- subset(soi, t > as.POSIXct("2000-01-01"))
##D oce.plot.ts(soirecent$t, soirecent$index)
## End(Not run)



cleanEx()
nameEx("read.lobo")
### * read.lobo

flush(stderr()); flush(stdout())

### Name: read.lobo
### Title: Read a LOBO File
### Aliases: read.lobo

### ** Examples

## Not run: 
##D library(oce)
##D uri <- paste("http://lobo.satlantic.com/cgi-bin/nph-data.cgi?",
##D   "min_date=20070220&max_date=20070305",
##D   "&x=date&",
##D   "y=current_across1,current_along1,nitrate,fluorescence,salinity,temperature&",
##D   "data_format=text",sep="")
##D lobo <- read.lobo(uri)
## End(Not run)



cleanEx()
nameEx("read.met")
### * read.met

flush(stderr()); flush(stdout())

### Name: read.met
### Title: Read a met File
### Aliases: read.met

### ** Examples

# The old csv format (no longer supplied by Environment Canada as of Oct 2019)
csv1 <- read.met(system.file("extdata", "test_met_vsn1.csv", package="oce"))

# The new "csv2" format (provided by Environment Canada as of Oct 2019)
csv2 <- read.met(system.file("extdata", "test_met_vsn2.csv", package="oce"))

# "xml2" format
if (requireNamespace("XML", quietly=TRUE))
    xml2 <- read.met(system.file("extdata", "test_met_xml2.xml", package="oce"))

## Not run: 
##D library(oce)
##D # Recreate data(met) and plot u(t) and v(t)
##D metFile <- download.met(id=6358, year=2003, month=9, destdir=".")
##D met <- read.met(metFile)
##D met <- oceSetData(met, "time", met[["time"]]+4*3600,
##D                   note="add 4h to local time to get UTC time")
##D plot(met)
## End(Not run)




cleanEx()
nameEx("read.oce")
### * read.oce

flush(stderr()); flush(stdout())

### Name: read.oce
### Title: Read an Oceanographic Data File
### Aliases: read.oce

### ** Examples

library(oce)
x <- read.oce(system.file("extdata", "ctd.cnv", package="oce"))
plot(x) # summary with TS and profiles
plotTS(x) # just the TS



cleanEx()
nameEx("read.odf")
### * read.odf

flush(stderr()); flush(stdout())

### Name: read.odf
### Title: Read an ODF file
### Aliases: read.odf

### ** Examples

library(oce)
#
# 1. Read a CTD cast made on the Scotian Shelf. Note that the file's metadata
# states that conductivity is in S/m, but it is really conductivity ratio,
# so we must alter the unit before converting to a CTD object. Note that
# read.odf() on this data file produces a warning suggesting that the user
# repair the unit, using the method outlined here.
odf <- read.odf(system.file("extdata", "CTD_BCD2014666_008_1_DN.ODF.gz", package="oce"))
ctd <- as.ctd(odf) ## so we can e.g. extract potential temperature
ctd[["conductivityUnit"]] <- list(unit=expression(), scale="")
#
# 2. Make a CTD, and plot (with span to show NS)
plot(ctd, span=500)
#
# 3. Highlight bad data on TS diagram. (Note that the eos
# is specified, because we will extract practical-salinity and
# UNESCO-defined potential temperatures for the added points.)
plotTS(ctd, type="o", eos="unesco") # use a line to show loops
bad <- ctd[["QCFlag"]]!=0
points(ctd[['salinity']][bad],ctd[['theta']][bad],col='red',pch=20)




cleanEx()
nameEx("read.topo")
### * read.topo

flush(stderr()); flush(stdout())

### Name: read.topo
### Title: Read a Topo File
### Aliases: read.topo

### ** Examples

## Not run: 
##D library(oce)
##D topoMaritimes <- read.topo("topoMaritimes.asc")
##D plot(topographyMaritimes)
## End(Not run)




cleanEx()
nameEx("read.woa")
### * read.woa

flush(stderr()); flush(stdout())

### Name: read.woa
### Title: Read a World Ocean Atlas NetCDF File
### Aliases: read.woa

### ** Examples

## Not run: 
##D ## Mean SST at 5-degree spatial resolution
##D tmn <- read.woa("/data/woa13/woa13_decav_t00_5dv2.nc", "t_mn")
##D imagep(tmn$longitude, tmn$latitude, tmn$t_mn[,,1], zlab="SST")
## End(Not run)



cleanEx()
nameEx("read.xbt")
### * read.xbt

flush(stderr()); flush(stdout())

### Name: read.xbt
### Title: Read an xbt file
### Aliases: read.xbt

### ** Examples

library(oce)
xbt <- read.oce(system.file("extdata", "xbt.edf", package="oce"))
summary(xbt)
plot(xbt)




cleanEx()
nameEx("read.xbt.edf")
### * read.xbt.edf

flush(stderr()); flush(stdout())

### Name: read.xbt.edf
### Title: Read a Sippican '.edf' format xbt file
### Aliases: read.xbt.edf

### ** Examples

library(oce)
xbt <- read.oce(system.file("extdata", "xbt.edf", package="oce"))
summary(xbt)
plot(xbt)




cleanEx()
nameEx("rescale")
### * rescale

flush(stderr()); flush(stdout())

### Name: rescale
### Title: Rescale values to lie in a given range
### Aliases: rescale

### ** Examples

library(oce)
# Fake tow-yow data
t <- seq(0, 600, 5)
x <- 0.5 * t
z <- 50 * (-1 + sin(2 * pi * t / 360))
T <- 5 + 10 * exp(z / 100)
palette <- oce.colorsViridis(100)
zlim <- range(T)
drawPalette(zlim=zlim, col=palette)
plot(x, z, type='p', pch=20, cex=3,
     col=palette[rescale(T, xlow=zlim[1], xhigh=zlim[2], rlow=1, rhigh=100)])



cleanEx()
nameEx("retime")
### * retime

flush(stderr()); flush(stdout())

### Name: retime
### Title: Adjust the time within Oce object
### Aliases: retime

### ** Examples

library(oce)
data(adv)
adv2 <- retime(adv,0,1e-4,as.POSIXct("2008-07-01 00:00:00",tz="UTC"))
plot(adv[["time"]], adv2[["time"]]-adv[["time"]], type='l')



cleanEx()
nameEx("rotateAboutZ")
### * rotateAboutZ

flush(stderr()); flush(stdout())

### Name: rotateAboutZ
### Title: Rotate velocity components within an oce object
### Aliases: rotateAboutZ

### ** Examples

library(oce)
par(mfcol=c(2, 3))
# adp (acoustic Doppler profiler)
data(adp)
plot(adp, which="uv")
mtext("adp", side=3, line=0, adj=1, cex=0.7)
adpRotated <- rotateAboutZ(adp, 30)
plot(adpRotated, which="uv")
mtext("adp rotated 30 deg", side=3, line=0, adj=1, cex=0.7)
# adv (acoustic Doppler velocimeter)
data(adv)
plot(adv, which="uv")
mtext("adv", side=3, line=0, adj=1, cex=0.7)
advRotated <- rotateAboutZ(adv, 125)
plot(advRotated, which="uv")
mtext("adv rotated 125 deg", side=3, line=0, adj=1, cex=0.7)
# cm (current meter)
data(cm)
plot(cm, which="uv")
mtext("cm", side=3, line=0, adj=1, cex=0.7)
cmRotated <- rotateAboutZ(cm, 30)
plot(cmRotated, which="uv")
mtext("cm rotated 30 deg", side=3, line=0, adj=1, cex=0.7)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("rsk")
### * rsk

flush(stderr()); flush(stdout())

### Name: rsk
### Title: Sample Rsk Dataset
### Aliases: rsk

### ** Examples

library(oce)
data(rsk)
## The object doesn't "know" it is CTD until told so
plot(rsk)
plot(as.ctd(rsk))




cleanEx()
nameEx("rskPatm")
### * rskPatm

flush(stderr()); flush(stdout())

### Name: rskPatm
### Title: Estimate Atmospheric Pressure in Rsk Object
### Aliases: rskPatm

### ** Examples

library(oce)
data(rsk)
print(rskPatm(rsk))




cleanEx()
nameEx("rskToc")
### * rskToc

flush(stderr()); flush(stdout())

### Name: rskToc
### Title: Decode table-of-contents File from a Rsk File
### Aliases: rskToc

### ** Examples

## Not run: 
##D table <- rskToc("/data/archive/sleiwex/2008/moorings/m05/adv/sontek_202h/raw",
##D from=as.POSIXct("2008-07-01 00:00:00", tz="UTC"),
##D     to=as.POSIXct("2008-07-01 12:00:00", tz="UTC"))
##D print(table)
## End(Not run)




cleanEx()
nameEx("runlm")
### * runlm

flush(stderr()); flush(stdout())

### Name: runlm
### Title: Calculate running linear models
### Aliases: runlm

### ** Examples


library(oce)

# Case 1: smooth a noisy signal
x <- 1:100
y <- 1 + x/100 + sin(x/5)
yn <- y + rnorm(100, sd=0.1)
L <- 4
calc <- runlm(x, y, L=L)
plot(x, y, type='l', lwd=7, col='gray')
points(x, yn, pch=20, col='blue')
lines(x, calc$y, lwd=2, col='red')

# Case 2: square of buoyancy frequency
data(ctd)
par(mfrow=c(1,1))
plot(ctd, which="N2")
rho <- swRho(ctd)
z <- swZ(ctd)
zz <- seq(min(z), max(z), 0.1)
N2 <- -9.8 / mean(rho) * runlm(z, rho, zz, deriv=1)
lines(N2, -zz, col='red')
legend("bottomright", lwd=2, bg="white",
       col=c("black", "red"),
       legend=c("swN2()", "using runlm()"))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("sealevelTuktoyaktuk")
### * sealevelTuktoyaktuk

flush(stderr()); flush(stdout())

### Name: sealevelTuktoyaktuk
### Title: Sea-level data set acquired in 1975 at Tuktoyaktuk
### Aliases: sealevelTuktoyaktuk

### ** Examples





cleanEx()
nameEx("secondsToCtime")
### * secondsToCtime

flush(stderr()); flush(stdout())

### Name: secondsToCtime
### Title: Time interval as colon-separated string
### Aliases: secondsToCtime

### ** Examples

library(oce)
cat("   10 s = ", secondsToCtime(10), "\n", sep="")
cat("   61 s = ", secondsToCtime(61), "\n", sep="")
cat("86400 s = ", secondsToCtime(86400), "\n", sep="")



cleanEx()
nameEx("section-class")
### * section-class

flush(stderr()); flush(stdout())

### Name: section-class
### Title: Class to Store Hydrographic Section Data
### Aliases: section-class

### ** Examples

library(oce)
data(section)
plot(section[['station', 1]])
pairs(cbind(z=-section[["pressure"]],T=section[["temperature"]],S=section[["salinity"]]))
## T profiles for first few stations in section, at common scale
par(mfrow=c(3,3))
Tlim <- range(section[["temperature"]])
ylim <- rev(range(section[["pressure"]]))
for (stn in section[["station", 1:9]])
    plotProfile(stn, xtype='potential temperature', ylim=ylim, Tlim=Tlim)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("section")
### * section

flush(stderr()); flush(stdout())

### Name: section
### Title: Hydrographic section
### Aliases: section

### ** Examples

## Not run: 
##D library(oce)
##D # Gulf Stream
##D data(section)
##D GS <- subset(section, 109<=stationId&stationId<=129)
##D GSg <- sectionGrid(GS, p=seq(0, 5000, 100))
##D plot(GSg, map.xlim=c(-80,-60))
## End(Not run)




cleanEx()
nameEx("sectionAddStation")
### * sectionAddStation

flush(stderr()); flush(stdout())

### Name: sectionAddStation
### Title: Add a CTD Station to a Section
### Aliases: sectionAddStation sectionAddCtd

### ** Examples

library(oce)
data(ctd)
ctd2 <- ctd
ctd2[["temperature"]] <- ctd2[["temperature"]] + 0.5
ctd2[["latitude"]] <- ctd2[["latitude"]] + 0.1
section <- as.section(c("ctd", "ctd2"))
ctd3 <- ctd
ctd3[["temperature"]] <- ctd[["temperature"]] + 1
ctd3[["latitude"]] <- ctd[["latitude"]] + 0.1
ctd3[["station"]] <- "Stn 3"
sectionAddStation(section, ctd3)




cleanEx()
nameEx("sectionGrid")
### * sectionGrid

flush(stderr()); flush(stdout())

### Name: sectionGrid
### Title: Grid a Section in Pressure Space
### Aliases: sectionGrid

### ** Examples

# Gulf Stream
library(oce)
data(section)
GS <- subset(section, 109<=stationId&stationId<=129)
GSg <- sectionGrid(GS, p=seq(0, 5000, 100))
plot(GSg, map.xlim=c(-80,-60))
# Show effects of various depth schemes
par(mfrow=c(3, 1))
default <- sectionGrid(GS)
approxML <- sectionGrid(GS, method="approxML")
standardDepths5 <- sectionGrid(GS, p=standardDepths(5))
plot(default, which="temperature", ztype="image", ylim=c(200,0))
mtext("default sectionGrid()")
plot(approxML, which="temperature", ztype="image", ylim=c(200,0))
mtext("sectionGrid(..., method=\"approxML\")")
plot(standardDepths5, which="temperature", ztype="image", ylim=c(200,0))
mtext("sectionGrid(..., p=standardDepths(5))")




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("sectionSmooth")
### * sectionSmooth

flush(stderr()); flush(stdout())

### Name: sectionSmooth
### Title: Smooth a Section
### Aliases: sectionSmooth

### ** Examples

# Unsmoothed (Gulf Stream)
library(oce)
data(section)
gs <- subset(section, 115<=stationId&stationId<=125)
par(mfrow=c(2, 2))

plot(gs, which="temperature")
mtext("unsmoothed")

# Spline
gsg <- sectionGrid(gs, p=seq(0, 5000, 100))
gsSpline <- sectionSmooth(gsg, "spline")
plot(gsSpline, which="temperature")
mtext("spline-smoothed")

# Barnes
gsBarnes <- sectionSmooth(gs, "barnes", xr=50, yr=200)
plot(gsBarnes, which="temperature")
mtext("Barnes-smoothed")

# Kriging
if (requireNamespace("automap",quietly=TRUE)&&requireNamespace("sp",quietly=TRUE)) {
 krig <- function(x, y, F, xg, xr, yg, yr) {
   xy <- data.frame(x=x/xr, y=y/yr)
   K <- automap::autoKrige(F~1, remove_duplicates=TRUE,
                           input_data=sp::SpatialPointsDataFrame(xy, data.frame(F)),
                           new_data=sp::SpatialPoints(expand.grid(xg/xr, yg/yr)))
   matrix(K$krige_output@data$var1.pred, nrow=length(xg), ncol=length(yg))
 }
 gsKrig <- sectionSmooth(gs, krig)
 plot(gsKrig, which="temperature")
 mtext("Kriging-smoothed")
}




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("sectionSort")
### * sectionSort

flush(stderr()); flush(stdout())

### Name: sectionSort
### Title: Sort a Section
### Aliases: sectionSort

### ** Examples

## Not run: 
##D # Eastern North Atlantic, showing Mediterranean water;
##D # sorting by longitude makes it easier to compare
##D # the map and the section.
##D library(oce)
##D data(section)
##D s <- sectionGrid(subset(section, -30 <= longitude))
##D ss <- sectionSort(ss, by="longitude")
##D plot(ss)
## End(Not run)




cleanEx()
nameEx("setFlags-adp-method")
### * setFlags-adp-method

flush(stderr()); flush(stdout())

### Name: setFlags,adp-method
### Title: Set data-quality flags within a adp object
### Aliases: setFlags,adp-method

### ** Examples

library(oce)
data(adp)

## Example 1: flag first 10 samples in a mid-depth bin of beam 1
i1 <- data.frame(1:20, 40, 1)
adpQC <- initializeFlags(adp, "v", 2)
adpQC <- setFlags(adpQC, "v", i1, 3)
adpClean1 <- handleFlags(adpQC, flags=list(3), actions=list("NA"))
par(mfrow=c(2, 1))
## Top: original, bottom: altered
plot(adp, which="u1")
plot(adpClean1, which="u1")

## Example 2: percent-good and error-beam scheme
v <- adp[["v"]]
i2 <- array(FALSE, dim=dim(v))
g <- adp[["g", "numeric"]]
# Thresholds on percent "goodness" and error "velocity"
G <- 25
V4 <- 0.45
for (k in 1:3)
    i2[,,k] <- ((g[,,k]+g[,,4]) < G) | (v[,,4] > V4)
adpQC2 <- initializeFlags(adp, "v", 2)
adpQC2 <- setFlags(adpQC2, "v", i2, 3)
adpClean2 <- handleFlags(adpQC2, flags=list(3), actions=list("NA"))
## Top: original, bottom: altered
plot(adp, which="u1")
plot(adpClean2, which="u1") # differs at 8h and 20h




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("setFlags-ctd-method")
### * setFlags-ctd-method

flush(stderr()); flush(stdout())

### Name: setFlags,ctd-method
### Title: Set data-quality flags within a ctd object
### Aliases: setFlags,ctd-method

### ** Examples

library(oce)
# Example 1: Range-check salinity
data(ctdRaw)
## Salinity and temperature range checks
qc <- ctdRaw
# Initialize flags to 2, meaning good data in the default
# scheme for handleFlags(ctd).
qc <- initializeFlags(qc, "salinity", 2)
qc <- initializeFlags(qc, "temperature", 2)
# Flag bad salinities as 4
oddS <- with(qc[["data"]], salinity < 25 | 40 < salinity)
qc <- setFlags(qc, name="salinity", i=oddS, value=4)
# Flag bad temperatures as 4
oddT <- with(qc[["data"]], temperature < -2 | 40 < temperature)
qc <- setFlags(qc, name="temperature", i=oddT, value=4)
# Compare results in TS space
par(mfrow=c(2, 1))
plotTS(ctdRaw)
plotTS(handleFlags(qc, flags=c(1, 3:9)))

# Example 2: Interactive flag assignment based on TS plot, using
# WHP scheme to define 'acceptable' and 'bad' codes
## Not run: 
##D options(eos="gsw")
##D data(ctd)
##D qc <- ctd
##D qc <- initializeFlagScheme(qc, "WHP CTD")
##D qc <- initializeFlags(qc, "salinity", 2)
##D Sspan <- diff(range(qc[["SA"]]))
##D Tspan <- diff(range(qc[["CT"]]))
##D n <- length(qc[["SA"]])
##D par(mfrow=c(1, 1))
##D plotTS(qc, type="o")
##D message("Click on bad points; quit by clicking to right of plot")
##D for (i in seq_len(n)) {
##D     xy <- locator(1)
##D     if (xy$x > par("usr")[2])
##D         break
##D     i <- which.min(abs(qc[["SA"]] - xy$x)/Sspan + abs(qc[["CT"]] - xy$y)/Tspan)
##D     qc <- setFlags(qc, "salinity", i=i, value=4)
##D     qc <- handleFlags(qc, flags=list(salinity=4))
##D     plotTS(qc, type="o")
##D }
## End(Not run)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("showMetadataItem")
### * showMetadataItem

flush(stderr()); flush(stdout())

### Name: showMetadataItem
### Title: Show metadata item
### Aliases: showMetadataItem

### ** Examples

library(oce)
data(ctd)
showMetadataItem(ctd, "ship", "ship")



cleanEx()
nameEx("siderealTime")
### * siderealTime

flush(stderr()); flush(stdout())

### Name: siderealTime
### Title: Convert a POSIXt time to a sidereal time
### Aliases: siderealTime

### ** Examples


t <- ISOdatetime(1978, 11, 13, 0, 0, 0, tz="UTC")
print(siderealTime(t))




cleanEx()
nameEx("snakeToCamel")
### * snakeToCamel

flush(stderr()); flush(stdout())

### Name: snakeToCamel
### Title: Convert each of a vector of strings from SNAKE_CASE to camelCase
### Aliases: snakeToCamel

### ** Examples

library(oce)
snakeToCamel("PARAMETER_DATA_MODE")   # "parameterDataMode"
snakeToCamel("PARAMETER")             # "parameter"
snakeToCamel("HISTORY_QCTEST")        # "historyQctest"
snakeToCamel("HISTORY_QCTEST", "QC")  # "historyQCTest"
snakeToCamel("PROFILE_DOXY_QC")       # "profileDoxyQc"
snakeToCamel("PROFILE_DOXY_QC", "QC") # "profileDoxyQC"



cleanEx()
nameEx("standardDepths")
### * standardDepths

flush(stderr()); flush(stdout())

### Name: standardDepths
### Title: Standard Oceanographic Depths
### Aliases: standardDepths

### ** Examples

depth  <- standardDepths()
depth1 <- standardDepths(1)
plot(depth, depth)
points(depth1, depth1, col=2, pch=20, cex=1/2)




cleanEx()
nameEx("sub-sub-adp-method")
### * sub-sub-adp-method

flush(stderr()); flush(stdout())

### Name: [[,adp-method
### Title: Extract Something from an adp Object
### Aliases: [[,adp-method

### ** Examples

data(adp)
# Tests for beam 1, distance bin 1, first 5 observation times
adp[["v"]][1:5,1,1]
adp[["a"]][1:5,1,1]
adp[["a", "numeric"]][1:5,1,1]
as.numeric(adp[["a"]][1:5,1,1]) # same as above




cleanEx()
nameEx("sub-sub-adv-method")
### * sub-sub-adv-method

flush(stderr()); flush(stdout())

### Name: [[,adv-method
### Title: Extract Something from an adv Object
### Aliases: [[,adv-method

### ** Examples

data(adv)
head(adv[["q"]])            # in raw form
head(adv[["q", "numeric"]]) # in numeric form




cleanEx()
nameEx("sub-sub-amsr-method")
### * sub-sub-amsr-method

flush(stderr()); flush(stdout())

### Name: [[,amsr-method
### Title: Extract Something From an amsr Object
### Aliases: [[,amsr-method

### ** Examples

# Histogram of SST values
library(oce)
data(amsr)
hist(amsr[["SST"]])




cleanEx()
nameEx("sub-sub-argo-method")
### * sub-sub-argo-method

flush(stderr()); flush(stdout())

### Name: [[,argo-method
### Title: Extract Something From an Argo Object
### Aliases: [[,argo-method

### ** Examples

data(argo)
# 1. show that dataset has 223 profiles, each with 56 levels
dim(argo[['temperature']])

# 2. show importance of focussing on data flagged 'good'
fivenum(argo[["salinity"]],na.rm=TRUE)
fivenum(argo[["salinity"]][argo[["salinityFlag"]]==1],na.rm=TRUE)




cleanEx()
nameEx("sub-sub-ctd-method")
### * sub-sub-ctd-method

flush(stderr()); flush(stdout())

### Name: [[,ctd-method
### Title: Extract Something From a ctd Object
### Aliases: [[,ctd-method

### ** Examples

data(ctd)
head(ctd[["temperature"]])




cleanEx()
nameEx("sub-sub-ladp-method")
### * sub-sub-ladp-method

flush(stderr()); flush(stdout())

### Name: [[,ladp-method
### Title: Extract Something From an ladp Object
### Aliases: [[,ladp-method

### ** Examples

data(ctd)
head(ctd[["temperature"]])




cleanEx()
nameEx("sub-sub-section-method")
### * sub-sub-section-method

flush(stderr()); flush(stdout())

### Name: [[,section-method
### Title: Extract Something From a Section Object
### Aliases: [[,section-method

### ** Examples

data(section)
length(section[["latitude"]])
length(section[["latitude", "byStation"]])
# Vector of all salinities, for all stations
Sv <- section[["salinity"]]
# List of salinities, grouped by station
Sl <- section[["salinity", "byStation"]]
# First station salinities
Sl[[1]]




cleanEx()
nameEx("sub-sub-topo-method")
### * sub-sub-topo-method

flush(stderr()); flush(stdout())

### Name: [[,topo-method
### Title: Extract Something From a Topo Object
### Aliases: [[,topo-method

### ** Examples

data(topoWorld)
dim(topoWorld[['z']])




cleanEx()
nameEx("sub-subset-ctd-method")
### * sub-subset-ctd-method

flush(stderr()); flush(stdout())

### Name: [[<-,ctd-method
### Title: Replace Parts of a ctd Object
### Aliases: [[<-,ctd-method

### ** Examples

data(ctd)
summary(ctd)
# Move the CTD profile a nautical mile north.
ctd[["latitude"]] <- 1/60 + ctd[["latitude"]] # acts in metadata
# Increase the salinity by 0.01.
ctd[["salinity"]] <- 0.01 + ctd[["salinity"]] # acts in data
summary(ctd)




cleanEx()
nameEx("sub-subset-section-method")
### * sub-subset-section-method

flush(stderr()); flush(stdout())

### Name: [[<-,section-method
### Title: Replace Parts of a Section Object
### Aliases: [[<-,section-method

### ** Examples

# 1. Change section ID from a03 to A03
data(section)
section[["sectionId"]]
section[["sectionId"]] <- toupper(section[["sectionId"]])
section[["sectionId"]]
# 2. Add a millidegree to temperatures at station 10
section[["station", 10]][["temperature"]] <-
    1e-3 + section[["station", 10]][["temperature"]]




cleanEx()
nameEx("subset-adp-method")
### * subset-adp-method

flush(stderr()); flush(stdout())

### Name: subset,adp-method
### Title: Subset an ADP Object
### Aliases: subset,adp-method

### ** Examples

library(oce)
data(adp)
# 1. Look at first part of time series, organized by time
earlyTime <- subset(adp, time < mean(range(adp[['time']])))
plot(earlyTime)

# 2. Look at first ten ensembles (AKA profiles)
en <- adp[["ensembleNumber"]]
firstTen <- subset(adp, ensembleNumber < en[11])
plot(firstTen)




cleanEx()
nameEx("subset-adv-method")
### * subset-adv-method

flush(stderr()); flush(stdout())

### Name: subset,adv-method
### Title: Subset an ADV Object
### Aliases: subset,adv-method

### ** Examples

library(oce)
data(adv)
plot(adv)
plot(subset(adv, time < mean(range(adv[['time']]))))




cleanEx()
nameEx("subset-amsr-method")
### * subset-amsr-method

flush(stderr()); flush(stdout())

### Name: subset,amsr-method
### Title: Subset an amsr Object
### Aliases: subset,amsr-method

### ** Examples

library(oce)
data(amsr) # see ?amsr for how to read and composite such objects
sub <- subset(amsr, -75 < longitude & longitude < -45)
sub <- subset(sub,   40 < latitude  &  latitude <  50)
plot(sub)
data(coastlineWorld)
lines(coastlineWorld[['longitude']], coastlineWorld[['latitude']])




cleanEx()
nameEx("subset-argo-method")
### * subset-argo-method

flush(stderr()); flush(stdout())

### Name: subset,argo-method
### Title: Subset an Argo Object
### Aliases: subset,argo-method subset.argo

### ** Examples

library(oce)
data(argo)

# Example 1: subset by time, longitude, and pressure
par(mfrow=c(2,2))
plot(argo)
plot(subset(argo, time > mean(time)))
plot(subset(argo, longitude > mean(longitude)))
plot(subset(argoGrid(argo), pressure > 500 & pressure < 1000), which=5)

# Example 2: restrict attention to delayed-mode profiles.
## Not run: 
##D par(mfrow=c(1, 1))
##D plot(subset(argo, dataMode == "D"))
## End(Not run)

# Example 3: contrast adjusted and unadjusted data
## Not run: 
##D par(mfrow=c(1, 2))
##D plotTS(argo)
##D plotTS(subset(argo, "adjusted"))
## End(Not run)

# Example 4. Subset by a polygon determined with locator()
## Not run: 
##D par(mfrow=c(1, 2))
##D plot(argo, which="map")
##D ## Can get a boundary with e.g. locator(4)
##D boundary <- list(x=c(-65, -40, -40, -65), y=c(65, 65, 45, 45))
##D argoSubset <- subset(argo, within=boundary)
##D plot(argoSubset, which="map")
## End(Not run)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("subset-cm-method")
### * subset-cm-method

flush(stderr()); flush(stdout())

### Name: subset,cm-method
### Title: Subset a CM Object
### Aliases: subset,cm-method

### ** Examples

library(oce)
data(cm)
plot(cm)
plot(subset(cm, time < mean(range(cm[['time']]))))




cleanEx()
nameEx("subset-coastline-method")
### * subset-coastline-method

flush(stderr()); flush(stdout())

### Name: subset,coastline-method
### Title: Subset a Coastline Object
### Aliases: subset,coastline-method subset.coastline

### ** Examples

library(oce)
data(coastlineWorld)
## Subset to a box centred on Nova Scotia, Canada
if (requireNamespace("sf")) {
    cl <- subset(coastlineWorld, -80<lon & lon<-50 & 30<lat & lat<60)
    ## The plot demonstrates that the trimming is as requested.
    plot(cl, clon=-65, clat=45, span=6000)
    rect(-80, 30, -50, 60, bg="transparent", border="red")
}



cleanEx()
nameEx("subset-ctd-method")
### * subset-ctd-method

flush(stderr()); flush(stdout())

### Name: subset,ctd-method
### Title: Subset a ctd Object
### Aliases: subset,ctd-method

### ** Examples

library(oce)
data(ctd)
plot(ctd)
## Example 1
plot(subset(ctd, pressure<10))
## Example 2
plot(subset(ctd, indices=1:10))




cleanEx()
nameEx("subset-echosounder-method")
### * subset-echosounder-method

flush(stderr()); flush(stdout())

### Name: subset,echosounder-method
### Title: Subset an Echosounder Object
### Aliases: subset,echosounder-method

### ** Examples

library(oce)
data(echosounder)
plot(echosounder)
plot(subset(echosounder, depth < 10))
plot(subset(echosounder, time < mean(range(echosounder[['time']]))))




cleanEx()
nameEx("subset-met-method")
### * subset-met-method

flush(stderr()); flush(stdout())

### Name: subset,met-method
### Title: Subset a met Object
### Aliases: subset,met-method

### ** Examples

library(oce)
data(met)
# Few days surrounding Hurricane Juan
plot(subset(met, time > as.POSIXct("2003-09-27", tz="UTC")))




cleanEx()
nameEx("subset-oce-method")
### * subset-oce-method

flush(stderr()); flush(stdout())

### Name: subset,oce-method
### Title: Subset an oce Object
### Aliases: subset,oce-method

### ** Examples

library(oce)
data(ctd)
# Select just the top 10 metres (pressure less than 10 dbar)
top10 <- subset(ctd, pressure < 10)
par(mfrow=c(1, 2))
plotProfile(ctd)
plotProfile(top10)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("subset-rsk-method")
### * subset-rsk-method

flush(stderr()); flush(stdout())

### Name: subset,rsk-method
### Title: Subset a Rsk Object
### Aliases: subset,rsk-method

### ** Examples

library(oce)
data(rsk)
plot(rsk)
plot(subset(rsk, time < mean(range(rsk[['time']]))))




cleanEx()
nameEx("subset-sealevel-method")
### * subset-sealevel-method

flush(stderr()); flush(stdout())

### Name: subset,sealevel-method
### Title: Subset a Sealevel Object
### Aliases: subset,sealevel-method

### ** Examples

library(oce)
data(sealevel)
plot(sealevel)
plot(subset(sealevel, time < mean(range(sealevel[['time']]))))




cleanEx()
nameEx("subset-section-method")
### * subset-section-method

flush(stderr()); flush(stdout())

### Name: subset,section-method
### Title: Subset a Section Object
### Aliases: subset,section-method subset.section

### ** Examples

library(oce)
data(section)

# Example 1. Stations within 500 km of the first station
starting <- subset(section, distance < 500)

# Example 2. Stations east of 50W
east <- subset(section, longitude > (-50))

# Example 3. Gulf Stream
GS <- subset(section, 109 <= stationId & stationId <= 129)

# Example 4. Only stations with more than 5 pressure levels
long <- subset(section, length(pressure) > 5)

# Example 5. Only stations that have some data in top 50 dbar
surfacing <- subset(section, min(pressure) < 50)

# Example 6. Similar to #4, but done in more detailed way
long <- subset(section,
   indices=unlist(lapply(section[["station"]],
                  function(s)
                    5 < length(s[["pressure"]]))))

# Example 7. Subset by a polygon determined with locator()
## Not run: 
##D par(mfrow=c(2, 1))
##D plot(section, which="map")
##D bdy <- locator(4) # choose a polygon near N. America
##D GS <- subset(section, within=bdy)
##D plot(GS, which="map")
## End(Not run)




cleanEx()
nameEx("subset-topo-method")
### * subset-topo-method

flush(stderr()); flush(stdout())

### Name: subset,topo-method
### Title: Subset a Topo Object
### Aliases: subset,topo-method

### ** Examples

## northern hemisphere
library(oce)
data(topoWorld)
plot(subset(topoWorld, latitude > 0))




cleanEx()
nameEx("subset-xbt-method")
### * subset-xbt-method

flush(stderr()); flush(stdout())

### Name: subset,xbt-method
### Title: Subset an xbt Object
### Aliases: subset,xbt-method

### ** Examples

library(oce)
data(xbt)
plot(xbt)
plot(subset(xbt, depth < mean(range(xbt[["depth"]]))))




cleanEx()
nameEx("summary-adv-method")
### * summary-adv-method

flush(stderr()); flush(stdout())

### Name: summary,adv-method
### Title: Summarize an ADV object
### Aliases: summary,adv-method

### ** Examples

library(oce)
data(adv)
summary(adv)




cleanEx()
nameEx("summary-argo-method")
### * summary-argo-method

flush(stderr()); flush(stdout())

### Name: summary,argo-method
### Title: Summarize an Argo Object
### Aliases: summary,argo-method summary.argo

### ** Examples

library(oce)
data(argo)
summary(argo)




cleanEx()
nameEx("summary-cm-method")
### * summary-cm-method

flush(stderr()); flush(stdout())

### Name: summary,cm-method
### Title: Summarize a CM Object
### Aliases: summary,cm-method

### ** Examples

library(oce)
data(cm)
summary(cm)




cleanEx()
nameEx("summary-ctd-method")
### * summary-ctd-method

flush(stderr()); flush(stdout())

### Name: summary,ctd-method
### Title: Summarize a ctd Object
### Aliases: summary,ctd-method

### ** Examples

library(oce)
data(ctd)
summary(ctd)




cleanEx()
nameEx("summary-lisst-method")
### * summary-lisst-method

flush(stderr()); flush(stdout())

### Name: summary,lisst-method
### Title: Summarize a LISST Object
### Aliases: summary,lisst-method

### ** Examples

library(oce)
data(lisst)
summary(lisst)




cleanEx()
nameEx("summary-lobo-method")
### * summary-lobo-method

flush(stderr()); flush(stdout())

### Name: summary,lobo-method
### Title: Summarize a LOBO Object
### Aliases: summary,lobo-method

### ** Examples


library(oce)
data(lobo)
summary(lobo)



cleanEx()
nameEx("summary-oce-method")
### * summary-oce-method

flush(stderr()); flush(stdout())

### Name: summary,oce-method
### Title: Summarize an oce Object
### Aliases: summary,oce-method

### ** Examples

o <- new("oce")
summary(o)



cleanEx()
nameEx("summary-rsk-method")
### * summary-rsk-method

flush(stderr()); flush(stdout())

### Name: summary,rsk-method
### Title: Summarize a Rsk Object
### Aliases: summary,rsk-method

### ** Examples

library(oce)
data(rsk)
summary(rsk)




cleanEx()
nameEx("summary-sealevel-method")
### * summary-sealevel-method

flush(stderr()); flush(stdout())

### Name: summary,sealevel-method
### Title: Summarize a Sealevel Object
### Aliases: summary,sealevel-method

### ** Examples

library(oce)
data(sealevel)
summary(sealevel)




cleanEx()
nameEx("summary-section-method")
### * summary-section-method

flush(stderr()); flush(stdout())

### Name: summary,section-method
### Title: Summarize a Section Object
### Aliases: summary,section-method summary.section

### ** Examples

library(oce)
data(section)
summary(section)




cleanEx()
nameEx("summary-tidem-method")
### * summary-tidem-method

flush(stderr()); flush(stdout())

### Name: summary,tidem-method
### Title: Summarize a Tidem Object
### Aliases: summary,tidem-method

### ** Examples

## Not run: 
##D library(oce)
##D data(sealevel)
##D tide <- tidem(sealevel)
##D summary(tide)
## End(Not run)




cleanEx()
nameEx("summary-topo-method")
### * summary-topo-method

flush(stderr()); flush(stdout())

### Name: summary,topo-method
### Title: Summarize A Topo Object
### Aliases: summary,topo-method

### ** Examples

library(oce)
data(topoWorld)
summary(topoWorld)




cleanEx()
nameEx("sunAngle")
### * sunAngle

flush(stderr()); flush(stdout())

### Name: sunAngle
### Title: Solar Angle as Function of Space and Time
### Aliases: sunAngle

### ** Examples


rise <- as.POSIXct("2011-03-03 06:49:00", tz="UTC") + 4*3600
set <- as.POSIXct("2011-03-03 18:04:00", tz="UTC") + 4*3600
mismatch <- function(lonlat)
{
    sunAngle(rise, lonlat[1], lonlat[2])$altitude^2 + sunAngle(set, lonlat[1], lonlat[2])$altitude^2
}
result <- optim(c(1,1), mismatch)
lon.hfx <- (-63.55274)
lat.hfx <- 44.65
dist <- geodDist(result$par[1], result$par[2], lon.hfx, lat.hfx)
cat(sprintf("Infer Halifax latitude %.2f and longitude %.2f; distance mismatch %.0f km",
            result$par[2], result$par[1], dist))




cleanEx()
nameEx("sunDeclinationRightAscension")
### * sunDeclinationRightAscension

flush(stderr()); flush(stdout())

### Name: sunDeclinationRightAscension
### Title: Sun Declination and Right Ascension
### Aliases: sunDeclinationRightAscension

### ** Examples

## Example 24.a in Meeus (1991) (page 158 PDF, 153 print)
time <- as.POSIXct("1992-10-13 00:00:00", tz="UTC")
a <- sunDeclinationRightAscension(time, apparent=TRUE)
expect_equal(a$declination, -7.78507,
             tol=            0.00004, scale=1)
expect_equal(a$rightAscension, -161.61919,
             tol=                 0.00003, scale=1)
b <- sunDeclinationRightAscension(time)
## check against previous results, to protect aginst code-drift errors
expect_equal(b$declination, -7.785464443,
             tol=            0.000000001, scale=1)
expect_equal(b$rightAscension, -161.6183305,
             tol=                 0.0000001, scale=1)




cleanEx()
nameEx("swAbsoluteSalinity")
### * swAbsoluteSalinity

flush(stderr()); flush(stdout())

### Name: swAbsoluteSalinity
### Title: Seawater absolute salinity, in GSW formulation
### Aliases: swAbsoluteSalinity

### ** Examples

## Not run: 
##D sa <- swAbsoluteSalinity(35.5, 300, 260, 16)
##D stopifnot(abs(35.671358392019094 - sa) < 00.000000000000010)
## End(Not run)




cleanEx()
nameEx("swAlphaOverBeta")
### * swAlphaOverBeta

flush(stderr()); flush(stdout())

### Name: swAlphaOverBeta
### Title: Ratio of seawater thermal expansion coefficient to haline
###   contraction coefficient
### Aliases: swAlphaOverBeta

### ** Examples

swAlphaOverBeta(40, 10, 4000, eos="unesco") # 0.3476




cleanEx()
nameEx("swCSTp")
### * swCSTp

flush(stderr()); flush(stdout())

### Name: swCSTp
### Title: Electrical conductivity ratio from salinity, temperature and
###   pressure
### Aliases: swCSTp

### ** Examples

expect_equal(1, swCSTp(35, T90fromT68(15), 0, eos="unesco")) # by definition of cond. ratio
expect_equal(1, swCSTp(34.25045, T90fromT68(15), 2000, eos="unesco"), tolerance=1e-7)
expect_equal(1, swCSTp(34.25045, T90fromT68(15), 2000, eos="gsw"), tolerance=1e-7)




cleanEx()
nameEx("swConservativeTemperature")
### * swConservativeTemperature

flush(stderr()); flush(stdout())

### Name: swConservativeTemperature
### Title: Seawater conservative temperature, in GSW formulation
### Aliases: swConservativeTemperature

### ** Examples

swConservativeTemperature(35,10,1000,188,4) # 9.86883




cleanEx()
nameEx("swDepth")
### * swDepth

flush(stderr()); flush(stdout())

### Name: swDepth
### Title: Water depth
### Aliases: swDepth

### ** Examples

d <- swDepth(10, 45)




cleanEx()
nameEx("swDynamicHeight")
### * swDynamicHeight

flush(stderr()); flush(stdout())

### Name: swDynamicHeight
### Title: Dynamic height of seawater profile
### Aliases: swDynamicHeight

### ** Examples

## Not run: 
##D library(oce)
##D data(section)
##D 
##D # Dynamic height and geostrophy
##D par(mfcol=c(2,2))
##D par(mar=c(4.5,4.5,2,1))
##D 
##D # Left-hand column: whole section
##D # (The smoothing lowers Gulf Stream speed greatly)
##D westToEast <- subset(section, 1<=stationId&stationId<=123)
##D dh <- swDynamicHeight(westToEast)
##D plot(dh$distance, dh$height, type='p', xlab="", ylab="dyn. height [m]")
##D ok <- !is.na(dh$height)
##D smu <- supsmu(dh$distance, dh$height)
##D lines(smu, col="blue")
##D f <- coriolis(section[["station", 1]][["latitude"]])
##D g <- gravity(section[["station", 1]][["latitude"]])
##D v <- diff(smu$y)/diff(smu$x) * g / f / 1e3 # 1e3 converts to m
##D plot(smu$x[-1], v, type='l', col="blue", xlab="distance [km]", ylab="velocity (m/s)")
##D 
##D # right-hand column: gulf stream region, unsmoothed
##D gs <- subset(section, 102<=stationId&stationId<=124)
##D dh.gs <- swDynamicHeight(gs)
##D plot(dh.gs$distance, dh.gs$height, type='b', xlab="", ylab="dyn. height [m]")
##D v <- diff(dh.gs$height)/diff(dh.gs$distance) * g / f / 1e3
##D plot(dh.gs$distance[-1], v, type='l', col="blue",
##D   xlab="distance [km]", ylab="velocity (m/s)")
## End(Not run)




cleanEx()
nameEx("swLapseRate")
### * swLapseRate

flush(stderr()); flush(stdout())

### Name: swLapseRate
### Title: Seawater lapse rate
### Aliases: swLapseRate

### ** Examples

lr <- swLapseRate(40, 40, 10000) # 3.255976e-4




cleanEx()
nameEx("swN2")
### * swN2

flush(stderr()); flush(stdout())

### Name: swN2
### Title: Squared buoyancy frequency for seawater
### Aliases: swN2

### ** Examples


library(oce)
data(ctd)
# Left panel: density
p <- ctd[["pressure"]]
ylim <- rev(range(p))
par(mfrow=c(1, 2), mar=c(3, 3, 1, 1), mgp=c(2, 0.7, 0))
plot(ctd[["sigmaTheta"]], p, ylim=ylim, type='l', xlab=expression(sigma[theta]))
# Right panel: N2, with default settings (black) and with df=2 (red)
N2 <- swN2(ctd)
plot(N2, p, ylim=ylim, xlab="N2 [1/s^2]", ylab="p", type="l")
lines(swN2(ctd, df=3), p, col=2)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("swPressure")
### * swPressure

flush(stderr()); flush(stdout())

### Name: swPressure
### Title: Water pressure
### Aliases: swPressure

### ** Examples

swPressure(9712.653, 30, eos="unesco") # 10000
swPressure(9712.653, 30, eos="gsw")    #  9998.863




cleanEx()
nameEx("swRho")
### * swRho

flush(stderr()); flush(stdout())

### Name: swRho
### Title: Seawater density
### Aliases: swRho

### ** Examples

library(oce)
# The numbers in the comments are the check values listed in reference 1;
# note that temperature in that reference was on the T68 scale, but that
# the present function works with the ITS-90 scale, so a conversion
# is required.
swRho(35, T90fromT68(5),      0, eos="unesco") # 1027.67547
swRho(35, T90fromT68(5),  10000, eos="unesco") # 1069.48914
swRho(35, T90fromT68(25),     0, eos="unesco") # 1023.34306
swRho(35, T90fromT68(25), 10000, eos="unesco") # 1062.53817




cleanEx()
nameEx("swRrho")
### * swRrho

flush(stderr()); flush(stdout())

### Name: swRrho
### Title: Density ratio
### Aliases: swRrho

### ** Examples

library(oce)
data(ctd)
u <- swRrho(ctd, eos="unesco")
g <- swRrho(ctd, eos="gsw")
p <- ctd[["p"]]
plot(u, p, ylim=rev(range(p)), type='l', xlab=expression(R[rho]))
lines(g, p, lty=2, col='red')
legend("topright", lty=1:2, legend=c("unesco", "gsw"), col=c("black", "red"))




cleanEx()
nameEx("swSCTp")
### * swSCTp

flush(stderr()); flush(stdout())

### Name: swSCTp
### Title: Practical salinity from electrical conductivity, temperature and
###   pressure
### Aliases: swSCTp

### ** Examples

# 1. Demonstrate agreement with test value in UNESCO documents
swSCTp(1, T90fromT68(15), 0, eos="unesco") # expect 35
# 2. Demonstrate agreement of gsw and unesco, S>2 case
swSCTp(1, T90fromT68(15), 0, eos="gsw") # again, expect 35
# 3. Demonstrate close values even in very brackish water
swSCTp(0.02, 10, 100, eos="gsw")    # 0.6013981
swSCTp(0.02, 10, 100, eos="unesco") # 0.6011721




cleanEx()
nameEx("swSTrho")
### * swSTrho

flush(stderr()); flush(stdout())

### Name: swSTrho
### Title: Seawater salinity from temperature and density
### Aliases: swSTrho

### ** Examples

swSTrho(10, 22, 0, eos="gsw") # 28.76285
swSTrho(10, 22, 0, eos="unesco") # 28.651625




cleanEx()
nameEx("swSigma")
### * swSigma

flush(stderr()); flush(stdout())

### Name: swSigma
### Title: Seawater density anomaly
### Aliases: swSigma

### ** Examples

library(oce)
swSigma(35, 13, 1000, longitude=300, latitude=30, eos="gsw") # 30.82374
swSigma(35, T90fromT68(13), 1000, eos="unesco")    # 30.8183




cleanEx()
nameEx("swSigmaT")
### * swSigmaT

flush(stderr()); flush(stdout())

### Name: swSigmaT
### Title: Seawater quasi-potential density anomaly
### Aliases: swSigmaT

### ** Examples

swSigmaT(35, 13, 1000, longitude=300, latitude=30, eos="gsw") # 26.39623
swSigmaT(35, T90fromT68(13), 1000, eos="unesco") # 26.39354




cleanEx()
nameEx("swSigmaTheta")
### * swSigmaTheta

flush(stderr()); flush(stdout())

### Name: swSigmaTheta
### Title: Seawater potential density anomaly
### Aliases: swSigmaTheta

### ** Examples

expect_equal(26.4212790994, swSigmaTheta(35, 13, 1000, eos="unesco"))




cleanEx()
nameEx("swSoundAbsorption")
### * swSoundAbsorption

flush(stderr()); flush(stdout())

### Name: swSoundAbsorption
### Title: Seawater sound absorption in dB/m
### Aliases: swSoundAbsorption

### ** Examples

## Fisher & Simmons (1977 table IV) gives 0.52 dB/km for 35 PSU, 5 degC, 500 atm
## (4990 dbar of water)a and 10 kHz
alpha <- swSoundAbsorption(35, 4, 4990, 10e3)

## reproduce part of Fig 8 of Francois and Garrison (1982 Fig 8)
f <- 1e3 * 10^(seq(-1,3,0.1)) # in KHz
plot(f/1000, 1e3*swSoundAbsorption(f, 35, 10, 0, formulation='fr'),
     xlab=" Freq [kHz]", ylab=" dB/km", type='l', log='xy')
lines(f/1000, 1e3*swSoundAbsorption(f, 0, 10, 0, formulation='fr'), lty='dashed')
legend("topleft", lty=c("solid", "dashed"), legend=c("S=35", "S=0"))




cleanEx()
nameEx("swSoundSpeed")
### * swSoundSpeed

flush(stderr()); flush(stdout())

### Name: swSoundSpeed
### Title: Seawater sound speed
### Aliases: swSoundSpeed

### ** Examples

swSoundSpeed(40, T90fromT68(40), 10000) # 1731.995 (p48 of Fofonoff + Millard 1983)




cleanEx()
nameEx("swSpecificHeat")
### * swSpecificHeat

flush(stderr()); flush(stdout())

### Name: swSpecificHeat
### Title: Seawater specific heat Source=
###   http://sam.ucsd.edu/sio210/propseawater/ppsw_fortran/ppsw.f check
###   value: cpsw = 3849.500 j/(kg deg. c) for s = 40 (ipss-78),
### Aliases: swSpecificHeat

### ** Examples

swSpecificHeat(40, T90fromT68(40), 10000, eos="unesco") # 3949.499




cleanEx()
nameEx("swSpice")
### * swSpice

flush(stderr()); flush(stdout())

### Name: swSpice
### Title: Seawater spiciness
### Aliases: swSpice

### ** Examples

## Contrast the two formulations.
library(oce)
data(ctd)
p <- ctd[["pressure"]]
plot(swSpice(ctd, eos="unesco"), p,
     xlim=c(-2.7, -1.5), ylim=rev(range(p)),
     xlab="Spice", ylab="Pressure (dbar)")
points(swSpice(ctd, eos="gsw"), p,col=2)
mtext("black=unesco, red=gsw")




cleanEx()
nameEx("swTFreeze")
### * swTFreeze

flush(stderr()); flush(stdout())

### Name: swTFreeze
### Title: Seawater freezing temperature
### Aliases: swTFreeze

### ** Examples

# 1. Test for a check-value given in reference 1. This value, -2.588567 degC,
# is in the 1968 temperature scale (IPTS-68), but swTFreeze reports
# in the newer ITS-90 scale, so we must convert before checking.
Tcheck <- -2.588567 # IPTS-68
T <- swTFreeze(salinity=40, pressure=500, eos="unesco")
expect_equal(Tcheck, T68fromT90(T), tolerance=1e-6)

# 2. Compare unesco and gsw formulations.
data(ctd)
p <- ctd[["pressure"]]
par(mfrow=c(1, 2), mar=c(3, 3, 1, 2), mgp=c(2, 0.7, 0))
plot(swTFreeze(ctd, eos="unesco"),
     p, xlab="unesco", ylim=rev(range(p)))
plot(swTFreeze(ctd, eos="unesco") - swTFreeze(ctd, eos="gsw"),
     p, xlab="unesco-gsw", ylim=rev(range(p)))




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("swTSrho")
### * swTSrho

flush(stderr()); flush(stdout())

### Name: swTSrho
### Title: Seawater temperature from salinity and density
### Aliases: swTSrho

### ** Examples

swTSrho(35, 23, 0, eos="unesco") # 26.11301




cleanEx()
nameEx("swThermalConductivity")
### * swThermalConductivity

flush(stderr()); flush(stdout())

### Name: swThermalConductivity
### Title: Seawater thermal conductivity
### Aliases: swThermalConductivity

### ** Examples

library(oce)
# Values in m^2/s, a unit that is often used instead of W/(m*degC).
swThermalConductivity(35, 10, 100) / (swRho(35,10,100) * swSpecificHeat(35,10,100)) # ocean
swThermalConductivity(0, 20, 0) / (swRho(0, 20, 0) * swSpecificHeat(0, 20, 0)) # lab
# Caldwell Table 1 gives 1478e-6 cal/(cm*sec*degC) at 31.5 o/oo, 10degC, 1kbar
joulePerCalorie <- 4.18400
cmPerM <- 100
swThermalConductivity(31.5,10,1000) / joulePerCalorie / cmPerM




cleanEx()
nameEx("swTheta")
### * swTheta

flush(stderr()); flush(stdout())

### Name: swTheta
### Title: Seawater potential temperature
### Aliases: swTheta

### ** Examples

library(oce)
## test value from Fofonoff et al., 1983
expect_equal(36.8818748026, swTheta(40, T90fromT68(40), 10000, 0, eos="unesco"))

# Example from a cross-Atlantic section
data(section)
stn <- section[['station', 70]]
plotProfile(stn, 'theta', ylim=c(6000, 1000))
lines(stn[['temperature']], stn[['pressure']], lty=2)
legend("bottomright", lty=1:2,
       legend=c("potential", "in-situ"),
       bg='white', title="Station 70")




cleanEx()
nameEx("swViscosity")
### * swViscosity

flush(stderr()); flush(stdout())

### Name: swViscosity
### Title: Seawater viscosity
### Aliases: swViscosity

### ** Examples

swViscosity(30, 10) # 0.001383779




cleanEx()
nameEx("test_met_csv1.csv")
### * test_met_csv1.csv

flush(stderr()); flush(stdout())

### Name: test_met_csv1.csv
### Title: Sample met dataset (CSV format prior to October 2019)
### Aliases: test_met_csv1.csv

### ** Examples

## Not run: 
##D read.met(system.file("extdata", "test_met_csv1.csv", package="oce"))
## End(Not run)




cleanEx()
nameEx("test_met_csv2.csv")
### * test_met_csv2.csv

flush(stderr()); flush(stdout())

### Name: test_met_csv2.csv
### Title: Sample met dataset (CSV format as of October 2019)
### Aliases: test_met_csv2.csv

### ** Examples

## Not run: 
##D read.met(system.file("extdata", "test_met_csv2.csv", package="oce"))
## End(Not run)




cleanEx()
nameEx("test_met_xml2.xml")
### * test_met_xml2.xml

flush(stderr()); flush(stdout())

### Name: test_met_xml2.xml
### Title: Sample met dataset (XML format as of October 2019)
### Aliases: test_met_xml2.xml

### ** Examples

## Not run: 
##D read.met(system.file("extdata", "test_met_xml2.xml", package="oce"))
## End(Not run)




cleanEx()
nameEx("threenum")
### * threenum

flush(stderr()); flush(stdout())

### Name: threenum
### Title: Calculate min, mean, and max values
### Aliases: threenum

### ** Examples

library(oce)
threenum(1:10)



cleanEx()
nameEx("tidalCurrent")
### * tidalCurrent

flush(stderr()); flush(stdout())

### Name: tidalCurrent
### Title: Tidal Current Dataset
### Aliases: tidalCurrent

### ** Examples

library(oce)
data(tidalCurrent)
par(mfrow=c(2, 1))
oce.plot.ts(tidalCurrent$time, tidalCurrent$u, ylab="u [m/s]")
abline(h=0, col=2)
oce.plot.ts(tidalCurrent$time, tidalCurrent$v, ylab="v [m/s]")
abline(h=0, col=2)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("tidem")
### * tidem

flush(stderr()); flush(stdout())

### Name: tidem
### Title: Fit a Tidal Model to a Timeseries
### Aliases: tidem

### ** Examples

library(oce)
# The demonstration time series from Foreman (1978),
# also used in T_TIDE (Pawlowicz, 2002).
data(sealevelTuktoyaktuk)
tide <- tidem(sealevelTuktoyaktuk)
summary(tide)

# AIC analysis
extractAIC(tide[["model"]])

# Fake data at M2
t <- seq(0, 10*86400, 3600)
eta <- sin(0.080511401 * t * 2 * pi / 3600)
sl <- as.sealevel(eta)
m <- tidem(sl)
summary(m)




cleanEx()
nameEx("tidemAstron")
### * tidemAstron

flush(stderr()); flush(stdout())

### Name: tidemAstron
### Title: Astronomical Calculations for Tidem
### Aliases: tidemAstron

### ** Examples

tidemAstron(as.POSIXct("2008-01-22 18:50:24"))




cleanEx()
nameEx("tidemVuf")
### * tidemVuf

flush(stderr()); flush(stdout())

### Name: tidemVuf
### Title: Nodal Modulation Calculations for Tidem
### Aliases: tidemVuf

### ** Examples

## Look up values for the M2 constituent in Halifax Harbour, Canada.
data("tidedata")
j  <- which(tidedata$const$name=="M2")
tidemVuf(t=as.POSIXct("2008-01-22 18:50:24"), j=j, lat=44.63)




cleanEx()
nameEx("timeToArgoJuld")
### * timeToArgoJuld

flush(stderr()); flush(stdout())

### Name: timeToArgoJuld
### Title: Convert time to Argo Julian Day (juld)
### Aliases: timeToArgoJuld

### ** Examples

timeToArgoJuld("2020-07-01")




cleanEx()
nameEx("topoInterpolate")
### * topoInterpolate

flush(stderr()); flush(stdout())

### Name: topoInterpolate
### Title: Interpolate Within a Topo Object
### Aliases: topoInterpolate

### ** Examples

library(oce)
data(topoWorld)
# "The Gully", approx. 400m deep, connects Gulf of St Lawrence with North Atlantic
topoInterpolate(45, -57, topoWorld)




cleanEx()
nameEx("topoWorld")
### * topoWorld

flush(stderr()); flush(stdout())

### Name: topoWorld
### Title: Global Topographic Dataset at Half-degree Resolution
### Aliases: topoWorld

### ** Examples

## Not run: 
##D library(oce)
##D data(topoWorld)
##D par(mfrow=c(2, 1))
##D plot(topoWorld, location=NULL)
##D imagep(topoWorld)
## End(Not run)




cleanEx()
nameEx("unabbreviateYear")
### * unabbreviateYear

flush(stderr()); flush(stdout())

### Name: unabbreviateYear
### Title: Determine year from various abbreviations
### Aliases: unabbreviateYear

### ** Examples

fullYear <- unabbreviateYear(c(99, 8, 108))



cleanEx()
nameEx("undriftTime")
### * undriftTime

flush(stderr()); flush(stdout())

### Name: undriftTime
### Title: Correct for drift in instrument clock
### Aliases: undriftTime

### ** Examples

## Not run: 
##D library(oce)
##D rbr011855 <- read.oce(
##D  "/data/archive/sleiwex/2008/moorings/m08/pt/rbr_011855/raw/pt_rbr_011855.dat")
##D d <- subset(rbr011855, time < as.POSIXct("2008-06-25 10:05:00"))
##D x <- undriftTime(d, 1)   # clock lost 1 second over whole experiment
##D summary(d)
##D summary(x)
## End(Not run)



cleanEx()
nameEx("unduplicateNames")
### * unduplicateNames

flush(stderr()); flush(stdout())

### Name: unduplicateNames
### Title: Rename duplicated character strings
### Aliases: unduplicateNames

### ** Examples

unduplicateNames(c("a", "b", "a", "c", "b"))
unduplicateNames(c("a", "b", "a", "c", "b"), style=2)



cleanEx()
nameEx("ungrid")
### * ungrid

flush(stderr()); flush(stdout())

### Name: ungrid
### Title: Extract (x, y, z) from (x, y, grid)
### Aliases: ungrid

### ** Examples

library(oce)
data(wind)
u <- interpBarnes(wind$x, wind$y, wind$z)
contour(u$xg, u$yg, u$zg)
U <- ungrid(u$xg, u$yg, u$zg)
points(U$x, U$y, col=oce.colorsViridis(100)[rescale(U$grid, rlow=1, rhigh=100)], pch=20)



cleanEx()
nameEx("unitFromString")
### * unitFromString

flush(stderr()); flush(stdout())

### Name: unitFromString
### Title: Decode units, from strings
### Aliases: unitFromString

### ** Examples

unitFromString("DB") # dbar



cleanEx()
nameEx("unwrapAngle")
### * unwrapAngle

flush(stderr()); flush(stdout())

### Name: unwrapAngle
### Title: Unwrap an angle that suffers modulo-360 problems
### Aliases: unwrapAngle

### ** Examples

library(oce)
true <- 355
a <- true + rnorm(100, sd=10)
a <- ifelse(a > 360, a - 360, a)
a2 <- unwrapAngle(a)
par(mar=c(3, 3, 5, 3))
hist(a, breaks=360)
abline(v=a2$mean, col="blue", lty="dashed")
abline(v=true, col="blue")
mtext("true (solid)\n estimate (dashed)", at=true, side=3, col="blue")
abline(v=mean(a), col="red")
mtext("mean", at=mean(a), side=3, col="red")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("utm2lonlat")
### * utm2lonlat

flush(stderr()); flush(stdout())

### Name: utm2lonlat
### Title: Convert UTM to Longitude and Latitude
### Aliases: utm2lonlat

### ** Examples

library(oce)
## Cape Split, in the Minas Basin of the Bay of Fundy
utm2lonlat(852863, 5029997, 19)




cleanEx()
nameEx("vectorShow")
### * vectorShow

flush(stderr()); flush(stdout())

### Name: vectorShow
### Title: Show some values from a vector
### Aliases: vectorShow

### ** Examples

vectorShow(pi)
vectorShow(volcano)
knot2mps <- 0.5144444
vectorShow(knot2mps, postscript="knots per m/s")
vectorShow("January", msg="The first month is")




cleanEx()
nameEx("velocityStatistics")
### * velocityStatistics

flush(stderr()); flush(stdout())

### Name: velocityStatistics
### Title: Report Statistics of adp or adv Velocities
### Aliases: velocityStatistics

### ** Examples

library(oce)
data(adp)
a <- velocityStatistics(adp)
print(a)
t <- seq(0, 2*pi, length.out=100)
theta <- a$ellipseAngle * pi / 180
y <- a$ellipseMajor * cos(t) * sin(theta) + a$ellipseMinor * sin(t) * cos(theta)
x <- a$ellipseMajor * cos(t) * cos(theta) - a$ellipseMinor * sin(t) * sin(theta)
plot(adp, which="uv+ellipse+arrow")
lines(x, y, col='blue', lty="dashed", lwd=5)
arrows(0, 0, a$uMean, a$vMean, lwd=5, length=1/10, col='blue', lty="dashed")




cleanEx()
nameEx("webtide")
### * webtide

flush(stderr()); flush(stdout())

### Name: webtide
### Title: Get a Tidal Prediction from a WebTide Database
### Aliases: webtide

### ** Examples

## Not run: 
##D ## needs WebTide at the system level
##D library(oce)
##D ## 1. prediction at Halifax NS
##D longitude <- -63.57
##D latitude <- 44.65
##D prediction <- webtide("predict", longitude=longitude, latitude=latitude)
##D mtext(sprintf("prediction at %fN %fE", latitude, longitude), line=0.75, side=3)
##D ## 2. map
##D webtide(lon=-63.57,lat=44.65,xlim=c(-64,-63),ylim=c(43.0,46))
## End(Not run)




cleanEx()
nameEx("window.oce")
### * window.oce

flush(stderr()); flush(stdout())

### Name: window.oce
### Title: Window an Oce Object by Time or Distance
### Aliases: window.oce

### ** Examples

library(oce)
data(adp)
plot(adp)
early <- window(adp, start="2008-06-26 00:00:00", end="2008-06-26 12:00:00")
plot(early)
bottom <- window(adp, start=0, end=20, which="distance")
plot(bottom)



cleanEx()
nameEx("write.ctd")
### * write.ctd

flush(stderr()); flush(stdout())

### Name: write.ctd
### Title: Write a CTD Data Object as a CSV File
### Aliases: write.ctd

### ** Examples

## Not run: 
##D library(oce)
##D data(ctd)
##D write.ctd(ctd, "ctd.csv")
##D d <- read.csv("ctd.csv")
##D plot(as.ctd(d$salinity, d$temperature, d$pressure))
## End(Not run)




cleanEx()
nameEx("xbt")
### * xbt

flush(stderr()); flush(stdout())

### Name: xbt
### Title: An XBT Object
### Aliases: xbt

### ** Examples

library(oce)
data(xbt)
summary(xbt)
plot(xbt)




cleanEx()
nameEx("xbt.edf")
### * xbt.edf

flush(stderr()); flush(stdout())

### Name: xbt.edf
### Title: Sample xbt dataset
### Aliases: xbt.edf

### ** Examples

## Not run: 
##D xbt <- read.oce(system.file("extdata", "xbt.edf", package="oce"))
## End(Not run)




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
