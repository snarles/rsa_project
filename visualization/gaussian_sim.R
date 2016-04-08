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

# source("visualization/source.R")
# cls_ind <- 1
# nrow(cls[[cls_ind]])
# smt <- get_surface(cls[[cls_ind]], cos_thres = 0.3, method = "randomForest")
# plot3d(cls[[cls_ind]], aspect = FALSE)
# points3d(smt, col = "red", size = 1)


xx <- cls[[cls_ind]]
mins <- apply(xx, 2, min) - 2
maxs <- apply(xx, 2, max) + 2
dat <- omask[mins[1]:maxs[1], mins[2]:maxs[2], mins[3]:maxs[3]]
# plot3d(which(dat, TRUE), aspect = FALSE)
# points3d(which(!dat, TRUE), col = "red", size = 1)
yvec <- as.factor(dat + 0)
yvec_num <- as.numeric(dat + 0)
xmat <- which(dat | TRUE, TRUE)
# yvec <- as.factor(dat[xmat] + 0)
colnames(xmat) <- NULL

#res <- randomForest(yvec ~ ., data = data.frame(yvec, xmat))
# res <- svm(yvec ~ ., data = data.frame(yvec, xmat), kernel = "radial", gamma = 0.9)
##res <- locfit(yvec ~ ., data = data.frame(yvec, xmat), family = "binomial", lfproc =)
# res <- loess(yvec_num ~ ., data = data.frame(yvec_num, xmat), span = 0.5)
# hist(predict(res, data = data.frame(yvec_num, xmat)))
# yvec2 <- (predict(res, data = data.frame(yvec_num, xmat)) > 0.2) + 0
# plot3d(xmat[yvec == 1, ], aspect = FALSE, size = 1)
# points3d(xmat[yvec2 == 1, ], col = "red", aspect = FALSE)
# table(yvec, yvec2)
prob_thres <- 0.66
pts_samp <-t(t(1.2 * pracma::rand(1e6, 3) - 0.1) * (maxs - mins))
# pts_samp <- rbind(xmat, pts_samp)
# plot3d(xmat, aspect = FALSE); points3d(pts_samp, col = "red", size = 1)
# y_samp <- (predict(res, data = data.frame(yvec = as.factor(rbinom(10000, 1, 0.5)), xmat = pts_samp)) > 0.2) + 0
y_samp <- knn(xmat, pts_samp, cl = yvec, k = 10, prob = TRUE)
y_prob <- attr(y_samp, "prob")
hist(y_prob)
min(y_prob)
# plot3d(pts_samp[y_samp == 1, ], size = 1 )
# plot3d(pts_samp[1:nrow(xmat), ][y_samp == 1, ])
plot3d(pts_samp[y_prob < prob_thres, ], size = 1, col = "red", aspect = FALSE)
points3d(xmat[yvec ==1, ])

plot3d(pts_samp[y_prob < prob_thres, ], size = 1, col = "red", aspect = FALSE)
