####
##  Simulate ROI shapes using gaussian random field
####

source("source1.R")
library(AlgDesign)
library(MASS)
library(rgl)
library(AnalyzeFMRI)
library(class)
library(locfit)
library(e1071)
library(randomForest)
gres <- 100
lbda <- 3
arr <- array(rnorm(gres^3), rep(gres, 3))
arr <- GaussSmoothArray(arr, ksize = 17, sigma = diag(rep(50, 3)))
omask <- arr > 0.04
## plot3d(which(arr > 0.04, TRUE))
cls <- get_cluster_inds2(omask, 20)
## delete any clusters which touch the edge
clsfilt <- sapply(cls, function(a) {
  mms <- c(apply(a, 2, min), apply(a, 2, max))
  !(1 %in% mms) && !(gres %in% mms)
})
cls <- cls[clsfilt]
plot3d(cls[[1]])

source("visualization/source.R")
cls_ind <- 1
xx <- cls[[cls_ind]]
mins <- apply(xx, 2, min) - 2
maxs <- apply(xx, 2, max) + 2
dat <- omask[mins[1]:maxs[1], mins[2]:maxs[2], mins[3]:maxs[3]]
smt <- get_surface2(dat, k=10, npoints = 1e5)
smt2 <- t(t(smt) + mins -2)
plot3d(xx, aspect = FALSE)
# plot3d(which(dat, TRUE), aspect = FALSE)
points3d(smt2, col = "red", size = 1)

