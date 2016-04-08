source("source1.R")

roi_dat <- readNIfTI('roi/zstat.nii.gz')
image(roi_dat)

## mask
roi_mask <- readNIfTI('roi/zstat_mask.nii.gz')
image(roi_mask)
mask_inds <- which(roi_mask == 1, arr.ind = TRUE)
mask_inds <- t(t(mask_inds) - c(27, 32, 23))
mask_inds <- mask_inds * 3
plot3d(mask_inds[mask_inds[, 1] < 0, ], aspect = FALSE)

xx <- mask_inds
plot3d(xx)
####
##  Filter down the xx matrix by eliminating low-magnitude points
####

mu <- colMeans(xx)
xx.c <- t(t(xx) - mu)
mags <- sqrt(rowSums(xx.c^2))
xx.n <- xx.c/mags
n <- nrow(xx)
filt <- rep(TRUE, n)
unvisited <- rep(TRUE, n)
cos_thres <- 0.001
while(sum(filt & unvisited) > 0) {
  ## determine next point to visit
  next_ind <- order(-mags + 100 * !(filt&unvisited))[1]
  mags[next_ind]
  unvisited[next_ind] <- FALSE
  ## find collection of close points to eliminate
  cds <- 1 - xx.n %*% xx.n[next_ind, ]
  elim_inds <- (cds < cos_thres) & unvisited
  filt[elim_inds] <- FALSE
  print(paste((1 - sum(filt & unvisited)/n) * 100 , "percent complete."))
}

xx2 <- xx.c[filt, ]
nrow(xx2)
plot3d(xx2, aspect = FALSE)

angles <- t(apply(xx2, 1, function(v) v/sqrt(sum(v^2)))); colnames(angles) <- NULL
mags2 <- mags[filt]
## interpolate mags as a function of angles

library(randomForest)
res <- randomForest(mags2 ~ ., data = data.frame(mags2 = mags2, a = angles))
angles_samp <- t(apply(pracma::randn(50000, 3), 1, function(v) v/sqrt(sum(v^2))))
mags_samp <- predict(res, data.frame(mags2 = 0, a = angles_samp))
xx_samp <- angles_samp * mags_samp
# plot3d(xx.c, col = "red", size = 2, aspect = FALSE)
# points3d(xx_samp)

plot3d(xx_samp, col = "black", size = 1, aspect = FALSE)
plot(xx_samp[, c(2, 1)], pch = ".", asp = 1, col = grey(0, 0.3))
plot(xx_samp[, c(1, 3)], pch = ".", asp = 1, col = grey(0, 0.3))
plot(xx_samp[, c(2, 3)], pch = ".", asp = 1, col = grey(0, 0.3))
plot(xx2[, c(2, 3)], pch = ".", asp = 1)

## convert back to original space

xx_decentered <- t(t(xx_samp) + mu)
plot3d(xx, size = 1, aspect = FALSE)
points3d(xx_decentered, col = "red")


## rescale

# mask_inds <- readRDS('visualization/brain_surface.rds')
mask_inds <- xx_decentered
# mask_inds <- t(t(mask_inds) - c(27, 32, 23))
# mask_inds <- mask_inds * 3
mask_inds <- mask_inds /2
mask_inds <- t(t(mask_inds) + c(46, 55, 46))

apply(mask_inds, 2, max)
apply(mask_inds, 2, min)

## save the sampled points!!

saveRDS(mask_inds, file = 'visualization/brain_surface.rds')
