####
##  Produce smoothed surfaces
####

library(randomForest)

get_surface <- function(xx, angles_samp = NULL, cos_thres = 0.1, method = c("randomForest", "loess"), ...) {
  method <- method[1]
  if (is.null(angles_samp)) {
    angles_samp <- t(apply(pracma::randn(5000, 3), 1, function(v) v/sqrt(sum(v^2))))
  }
  mu <- colMeans(xx)
  xx.c <- t(t(xx) - mu)
  mags <- sqrt(rowSums(xx.c^2))
  xx.n <- xx.c/mags
  n <- nrow(xx)
  filt <- rep(TRUE, n)
  unvisited <- rep(TRUE, n)
  while(sum(filt & unvisited) > 0) {
    ## determine next point to visit
    next_ind <- order(-mags + 100 * !(filt&unvisited))[1]
    mags[next_ind]
    unvisited[next_ind] <- FALSE
    ## find collection of close points to eliminate
    cds <- 1 - xx.n %*% xx.n[next_ind, ]
    elim_inds <- (cds < cos_thres) & unvisited
    filt[elim_inds] <- FALSE
    # print(paste((1 - sum(filt & unvisited)/n) * 100 , "percent complete."))
  }
  xx2 <- xx.c[filt, ]
  angles <- t(apply(xx2, 1, function(v) v/sqrt(sum(v^2)))); colnames(angles) <- NULL
  mags2 <- mags[filt]
  if (method == "randomForest") {
    res <- randomForest(mags2 ~ ., data = data.frame(mags2 = mags2, a = angles), ...)
  }
  if (method == "loess") {
    res <- loess(mags2 ~ ., data = data.frame(mags2 = mags2, a = angles), ...)
  }
  mags_samp <- predict(res, data.frame(mags2 = 0, a = angles_samp))
  xx_samp <- angles_samp * mags_samp
  xx_decentered <- t(t(xx_samp) + mu)
  xx_decentered
}

