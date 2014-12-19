##OLD ## http://r.789695.n4.nabble.com/question-re-error-message-package-error-quot-functionName-quot-not-resolved-from-current-namespace-td4663892.html
##OLD 
##OLD teos10 <- function(name, a1, a2, a3, a4, lib)
##OLD {
##OLD     if (missing(lib)) {
##OLD         lib <- paste(.libPaths(), "/teos10/libs/libgswteos10.so", sep="")
##OLD         ## FIXME: test for Windows machine and set up 'lib' appropriately 
##OLD     }
##OLD     if (missing(name))
##OLD         stop("a function name must be given, and it must be in lower case letters, e.g. \"gsw_sa_from_sp\"")
##OLD     if (missing(a1))
##OLD         stop("must provide a1")
##OLD     ## FIXME: later, can count when the missing list starts, and use that to figure out which C function to call
##OLD     if (missing(a2))
##OLD         stop("must provide a2")
##OLD     if (missing(a3)) {
##OLD         args <- 2
##OLD     } else {
##OLD         if (missing(a4)) {
##OLD             args <- 3 
##OLD         } else {
##OLD             args <- 4
##OLD         }                              # FIXME: this check on the number args is inelegant
##OLD     }
##OLD     dim <- dim(a1)
##OLD     a1 <- as.vector(a1)
##OLD     n <- length(a1)
##OLD     a2 <- as.vector(a2)
##OLD     ## FIXME: should use rep() to make lengths match, perhaps (surely this depends on teos calls though,
##OLD     ## e.g. consider pref)
##OLD     if (length(a2) != n)
##OLD         stop("length(a2) must match length(a1)")
##OLD     good <- is.finite(a1) & is.finite(a2)
##OLD     if (args > 2) {
##OLD         if (length(a3) != n)
##OLD             stop("length(a3) must match length(a1)")
##OLD         a3 <- as.vector(a3)
##OLD         good <- good & is.finite(a3)
##OLD     }
##OLD     if (args > 3) {
##OLD         if (length(a4) != n)
##OLD             stop("length(a4) must match length(a1)")
##OLD         a4 <- as.vector(a4)
##OLD         good <- good & is.finite(a4)
##OLD     }
##OLD     rval <- rep(NA, n)
##OLD     ngood <- sum(good)
##OLD     if (args == 2) {
##OLD         rval[good] <- .C("gsw2a", as.character(lib), as.character(name),
##OLD                          as.integer(ngood),
##OLD                          as.double(a1[good]),
##OLD                          as.double(a2[good]),
##OLD                          rval=double(ngood), NAOK=TRUE, PACKAGE="teos10")$rval
##OLD     } else if (args == 3) {
##OLD         rval[good] <- .C("gsw3a", as.character(lib), as.character(name),
##OLD                          as.integer(ngood),
##OLD                          as.double(a1[good]),
##OLD                          as.double(a2[good]),
##OLD                          as.double(a3[good]),
##OLD                          rval=double(ngood), NAOK=TRUE, PACKAGE="teos10")$rval
##OLD     } else if (args == 4) {
##OLD         rval[good] <- .C("gsw4a", as.character(lib), as.character(name),
##OLD                          as.integer(ngood),
##OLD                          as.double(a1[good]),
##OLD                          as.double(a2[good]),
##OLD                          as.double(a3[good]),
##OLD                          as.double(a4[good]),
##OLD                          rval=double(ngood), NAOK=TRUE, PACKAGE="teos10")$rval
##OLD     }
##OLD     rval[rval == 9e15] <- NA
##OLD     dim(rval) <- dim
##OLD     rval
##OLD }
##OLD 
