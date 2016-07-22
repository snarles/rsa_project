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

ind <- 7
mat <- cbind(paramz[[ind]][, -3], stypz[[ind]][, 1])
plot3d(mat)

