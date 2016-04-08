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
nrow(cls[[cls_ind]])
smt <- get_surface(cls[[cls_ind]], cos_thres = 0.3, method = "randomForest")
plot3d(cls[[cls_ind]], aspect = FALSE)
points3d(smt, col = "red", size = 1)


xx <- cls[[cls_ind]]
mins <- apply(xx, 2, min)
maxs <- apply(xx, 2, max)
dat <- omask[mins[1]:maxs[1], mins[2]:maxs[2], mins[3]:maxs[3]]
# plot3d(which(dat, TRUE), aspect = FALSE)
# points3d(which(!dat, TRUE), col = "red", size = 1)
yvec <- as.factor(dat + 0)
xmat <- which(dat | TRUE, TRUE)
# yvec <- as.factor(dat[xmat] + 0)
colnames(xmat) <- NULL

#res <- randomForest(yvec ~ ., data = data.frame(yvec, xmat))
res <- svm(yvec ~ ., data = data.frame(yvec, xmat), kernel = "radial", gamma = 0.9)
##res <- locfit(yvec ~ ., data = data.frame(yvec, xmat), family = "binomial", lfproc =)

yvec2 <- predict(res, data = data.frame(yvec, xmat))

plot3d(xmat[yvec == 1, ], aspect = FALSE, size = 1)
points3d(xmat[yvec2 == 1, ], col = "red", aspect = FALSE)

table(yvec, yvec2)

pts_samp <-t(t(pracma::rand(10000, 3)) * (maxs - mins)) + 1
pts_samp <- rbind(xmat, pts_samp)
# plot3d(xmat, aspect = FALSE); points3d(pts_samp, col = "red", size = 1)
y_samp <- predict(res, data = data.frame(yvec = as.factor(rbinom(10000, 1, 0.5)), xmat = pts_samp))
# y_samp <- knn(xmat, pts_samp, cl = yvec, k = 3)


plot3d(pts_samp[y_samp == 1, ])
plot3d(pts_samp[1:nrow(xmat), ][y_samp == 1, ])


y_samp <- predict(res, data = data.frame(yvec = as.factor(rbinom(10000, 1, 0.5)), xmat = pts_samp), type = "prob")
plot3d(pts_samp[y_samp[, 2] > 0.5, ])



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
