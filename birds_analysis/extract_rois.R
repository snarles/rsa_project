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


res0 <- oro.nifti::readNIfTI("DavisPoldrack_Birds_Archive/betaseries_MNI/bg_image.nii.gz")

rescovs <- list()

for (f in fnos) {
  fpath <- paste0("DavisPoldrack_Birds_Archive/betaseries_MNI/", f, "_lsone_MNI.nii.gz")
  res <- oro.nifti::readNIfTI(fpath)
  dim(res0) # [1]  91 109  91
  dim(res) # [1]  91 109  91 108
  
  gplot(res0[, , 30])
  gplot(res[, , 30, 100])
  
  plot(params$angle, res[35, 35, 30, ])
  plot(params$angle, params$height)
  
  x1 <- params$angle/sd(params$angle)
  x2 <- params$height/sd(params$height)
  
  resvar <- apply(res, 1:3, var)
  rescov1 <- apply(res, 1:3, function(v) cov(v, x1))
  rescov2 <- apply(res, 1:3, function(v) cov(v, x2))
  rescov <- sqrt(rescov1^2 + rescov2^2)
  rescovs[[f]] <- rescov
}
saveRDS(rescovs, file = "DavisPoldrack_Birds_Archive/rescovs.rds")

dim(resvar)
gplot(resvar[, , 30])





gplot(rescov1[, , 30])
gplot(rescov2[, , 30])


gplot(rescov[, , 30])

summary(as.numeric(rescov))
hist(as.numeric(rescov), breaks = 40)