####
##  analysis of gambling data
####

library(pracma)
library(lineId)
cls <- readRDS("roi/cl_inds.rds")
params <- readRDS("roi/params.rds")
dat <- readRDS("roi/data.rds")

## delete bad voxels
goodvox <- which(colSums(dat == 0) == 0)
cls <- cls[names(goodvox)[-(1:3)], ]
dat <- dat[, goodvox]

rdat <- dat[, -(1:3)]
hdat <- dat[, 1:3]
nrois <- max(cls[, "clus"])




####
##  Average multiple subjects
####

avg_dat <- 1/16 * repmat(eye(48), 1, 16) %*% rdat[, ]
dim(avg_dat) # [1]   48 5082
avg_dat <- scale(avg_dat)
params2 <- cbind(1, scale(params))
Xmat <- repmat(params2, 3, 1)
Ball <- fit_ridge_kernel(Xmat, avg_dat, lambda = 0.001)
Yhat <- Xmat %*% Ball      
dim(Yhat)
resid <- avg_dat - Yhat
resid_ses <- apply(resid, 2, sd)
hist(resid_ses)
colnames(Ball) <- colnames(avg_dat)
Bscale <- t(t(Ball/resid_ses))

####
##  Compute overall feature-distance matrix and roi-specific feature-distance matrices
####



