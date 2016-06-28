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
dim(res0) # [1]  91 109  91
dim(res) # [1]  91 109  91 108

gplot(res0[, , 30])
gplot(res[, , 30, 100])

plot(params$angle, res[35, 35, 30, ])
plot(params$angle, params$height)

x1 <- params$angle/sd(params$angle)
x2 <- params$height/sd(params$height)

resvar <- apply(res, 1:3, var)
dim(resvar)
gplot(resvar[, , 30])



rescov1 <- apply(res, 1:3, function(v) cov(v, x1))
rescov2 <- apply(res, 1:3, function(v) cov(v, x2))

gplot(rescov1[, , 30])
gplot(rescov2[, , 30])

rescov <- sqrt(rescov1^2 + rescov2^2)
gplot(rescov[, , 30])

summary(as.numeric(rescov))
hist(as.numeric(rescov), breaks = 40)

ids <- 10 * 1:9
for (z in ids)
  gplot(rescov[, , z])

for (z in ids)
  gplot(rescov[, , z] > 10)

for (z in ids)
  gplot(rescov[, , z] > 20)


for (z in ids)
  gplot(rescov[, , z] > 30)



