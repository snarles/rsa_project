####
##  Simulate ROI shapes using gaussian random field
####

source("source1.R")
library(AlgDesign)
library(MASS)
library(rgl)
library(AnalyzeFMRI)
gres <- 100
lbda <- 3
arr <- array(rnorm(gres^3), rep(gres, 3))
arr <- GaussSmoothArray(arr, ksize = 17, sigma = diag(rep(50, 3)))
## plot3d(which(arr > 0.04, TRUE))
cls <- get_cluster_inds2(arr > 0.04, 20)
## delete any clusters which touch the edge
clsfilt <- sapply(cls, function(a) {
  mms <- c(apply(a, 2, min), apply(a, 2, max))
  !(1 %in% mms) && !(gres %in% mms)
})
cls <- cls[clsfilt]
plot3d(cls[[1]])

####
##  Eliminate interior points
####

cosineDist <- function(x){
  1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
}

xx <- cls[[1]]
mu <- colMeans(xx)
xx.c <- t(t(xx) - mu)
dm <- cosineDist(xx.c)
mags <- sqrt(rowSums(xx.c^2))

## for every point, check if there is a larger point within some threshold
thres <- 0.01
filt <- logical(nrow(xx))
for (i in 1:nrow(xx)) {
  nbrs <- which(dm[i, ] < thres)
  filt[i] <- (max(mags[nbrs]) <= mags[i])
}
# plot3d(xx)
# points3d(xx[filt, ], col = "red", size = 5)


# plot3d(xx[filt, ])
# for (i in 1:sum(filt)) {
#   lines3d(rbind(xx[filt, ][i, ], mu))
# }


xx2 <- xx.c[filt, ]
angles <- t(apply(xx2, 1, function(v) v/sqrt(sum(v^2)))); colnames(angles) <- NULL
mags2 <- mags[filt]
## interpolate mags as a function of angles
angles_samp <- t(apply(pracma::randn(1000, 3), 1, function(v) v/sqrt(sum(v^2))))
res <- loess(mags2 ~ ., data = data.frame(mags2 = mags2, a = angles), span = 0.3)
mags_samp <- predict(res, data.frame(mags2 = 0, a = angles_samp))
xx_samp <- angles_samp * mags_samp
plot3d(xx.c, col = "red", size = 2)
points3d(xx_samp)
