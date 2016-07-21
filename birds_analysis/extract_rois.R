library(oro.nifti)
library(AnalyzeFMRI)
library(fmri)
library(pracma)

gplot <- function(sl) {
  image(t(flipud(sl)), col = gray.colors(10))
}

params <- read.table("DavisPoldrack_Birds_Archive/attribute_files/02_angleheight.txt", header = FALSE, sep = " ")
colnames(params) <- c("angle", "height", "run")


fnos <- c("002", "003", "004", "005", "006", "007", "008", "009", "011",
          "012", "013", "014", "017")

fnos2 <- c("02", "03", "04", "05", "06", "07", "08", "09", "11",
           "12", "13", "14", "17")

res0 <- oro.nifti::readNIfTI("DavisPoldrack_Birds_Archive/betaseries_MNI/bg_image.nii.gz")

# rescovs <- readRDS("DavisPoldrack_Birds_Archive/rescovs.rds")

rescovs <- list()
paramz <- list()

for (f in fnos2) {
  fpath2 <- paste0("DavisPoldrack_Birds_Archive/attribute_files/", f, "_angleheight.txt")
  params <- read.table(fpath2, header = FALSE, sep = " ")
  colnames(params) <- c("angle", "height", "run")
  paramz[[f]] <- params
}

for (f in fnos2) {
  fpath <- paste0("DavisPoldrack_Birds_Archive/betaseries_MNI/0", f, "_lsone_MNI.nii.gz")
  res <- oro.nifti::readNIfTI(fpath)
  fpath2 <- paste0("DavisPoldrack_Birds_Archive/attribute_files/", f, "_angleheight.txt")
  params <- read.table(fpath2, header = FALSE, sep = " ")
  colnames(params) <- c("angle", "height", "run")
  x1 <- params$angle/sd(params$angle)
  x2 <- params$height/sd(params$height)

  resvar <- apply(res, 1:3, var)
  rescov1 <- apply(res, 1:3, function(v) cov(v, x1))
  rescov2 <- apply(res, 1:3, function(v) cov(v, x2))
  rescov <- sqrt(rescov1^2 + rescov2^2)
  rescovs[[f]] <- rescov
}
saveRDS(rescovs, file = "DavisPoldrack_Birds_Archive/rescovs.rds")

length(rescovs)
names(rescovs)

hist(rescovs[[1]])


gplot(rescovs[[1]][, , 70] > 25)
gplot(rescovs[[2]][, , 70] > 25)
gplot(rescovs[[3]][, , 70] > 25)
gplot(rescovs[[4]][, , 70] > 25)
gplot(rescovs[[5]][, , 70] > 25)
gplot(rescovs[[6]][, , 70] > 25)
gplot(rescovs[[7]][, , 70] > 25)
gplot(rescovs[[8]][, , 70] > 25)
gplot(rescovs[[9]][, , 70] > 25)
gplot(rescovs[[10]][, , 70] > 25)
gplot(rescovs[[11]][, , 70] > 25)
gplot(rescovs[[12]][, , 70] > 25)
gplot(rescovs[[13]][, , 70] > 25)


