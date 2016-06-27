library(oro.nifti)
library(AnalyzeFMRI)
library(fmri)
library(pracma)

gplot <- function(sl) {
  image(t(flipud(sl)), col = gray.colors(10))
}

params <- read.table("DavisPoldrack_Birds_Archive/attribute_files/02_angleheight.txt", header = FALSE, sep = " ")
colnames(params) <- c("angle", "height", "run")


res0 <- oro.nifti::readNIfTI("DavisPoldrack_Birds_Archive/betaseries_MNI/bg_image.nii.gz")
res <- oro.nifti::readNIfTI("DavisPoldrack_Birds_Archive/betaseries_MNI/002_lsone_MNI.nii.gz")
dim(res0)
dim(res) # [1]  91 109  91 108

gplot(res0[, , 30])
gplot(res[, , 30, 100])

plot(params$angle, res[35, 35, 30, ])
plot(params$angle, params$height)

resvar <- apply(res, 1:3, var)
dim(resvar)
gplot(resvar[, , 30])
