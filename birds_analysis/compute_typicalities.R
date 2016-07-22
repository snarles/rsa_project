library(rgl)

fnos <- c("002", "003", "004", "005", "006", "007", "008", "009", "011",
          "012", "013", "014", "017")

fnos2 <- c("02", "03", "04", "05", "06", "07", "08", "09", "11",
           "12", "13", "14", "17")

stypz <- list()
paramz <- list()

for (f in fnos2) {
  fpath2 <- paste0("DavisPoldrack_Birds_Archive/attribute_files/", f, "_angleheight.txt")
  params <- read.table(fpath2, header = FALSE, sep = " ")
  colnames(params) <- c("angle", "height", "run")
  paramz[[f]] <- params
}

for (f in fnos) {
  fpath <- paste0("DavisPoldrack_Birds_Archive/attribute_files/", f, "_subject_typ.txt")
  styps <- read.table(fpath, header = FALSE, sep = " ")
  stypz[[f]] <- styps
}

for (ind in 1:length(fnos)) {
  pdf(paste0("birds_analysis/pplots/param", ind, ".pdf"))
  plot(paramz[[ind]][, -3])
  title(ind)
  dev.off()
}

ind <- 1
mat <- cbind(paramz[[ind]][, -3], stypz[[ind]][, 1])
plot3d(mat)


## OBSERVATIONS
##  Ind Sub  Params SubTyp  
##    1   2       A
##    2   3       B
##    3   4       A
##    4   5       C
##    5   6       B
##    6   7       D
##    7   8       C
##    8   9       B
##    9   11      C
##   10   12      A
##   11   13      B
##   12   14      C
##   13   17      B
