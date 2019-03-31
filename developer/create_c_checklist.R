## Extract function names from the C source. The method is brittle relative to
## changes in the C source format.

## We do a step-by-step dismantingly of the strings, because might make
## it easier to change the method, if the C function format changes.
source <- "~/git/GSW-C/gswteos-10.h"
C <- readLines(source)
C2 <- C[grep("extern", C)]
C3 <- C2[grep('extern "C" \\{', C2, invert=TRUE)]
C4 <- gsub("^[ ]*extern[ ]*((void)|(int)|(double))*[ ]*", "", C3)
C5 <- gsub("^\\*", "", C4)
C6 <- gsub("\\(.*$", "", C5)
C7 <- sort(C6)
cat("C functions in '", source, "' (", length(C7), ") are:\n* [ ] ",
        paste(C7, collapse="\n* [ ] "), "\n")
               
