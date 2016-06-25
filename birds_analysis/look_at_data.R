library(oro.nifti)
library(AnalyzeFMRI)
library(fmri)
library(pracma)

gplot <- function(sl) {
  image(t(flipud(sl)), col = gray.colors(10))
}

res <- oro.nifti::readNIfTI("DavisPoldrack_Birds_Archive/betaseries_MNI/bg_image.nii.gz")
res <- oro.nifti::readNIfTI("DavisPoldrack_Birds_Archive/betaseries_MNI/002_lsone_MNI.nii.gz")
dim(res)

gplot(res[, , 30, 100])
