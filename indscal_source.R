####
##  INDSCAL analysis ("Multiple MDS")
####

library(smacof)

sqrt2 <- function(v) sqrt(pmax(v, 0))

#' @param dists Squared euclidean distances
indscal_routine <- function(dists, p = 2) {
  n <- dim(dists[[1]])[1]
  k <- length(dists)
  of <- function(x, ds, dists) {
    nms <- lapply(1:k, function(i) colSums(t(x^2) * ds[i, ]))
    ips <- lapply(1:k, function(i) x %*% diag(ds[i, ]) %*% t(x))
    dhats <- lapply(1:k, function(i) {
      v <- t(nms[[i]])
      repmat(v, n, 1) + repmat(t(v), 1, n) - 2 * ips[[i]]
    })
    diffs <- sapply(1:k, function(i) f2(dists[[i]], dhats[[i]]))
    sum(diffs)
  }
  
  of_conv <- function(xx) {
    xx <- matrix(xx, n + k, p)
    x <- xx[1:n, ]
    ds <- xx[-(1:n), ]
    of(x, ds, dists)
  }
  
  dists2 <- lapply(dists, sqrt2)
  res <- indscal(dists2, ndim = p , type="ratio", verbose=FALSE)
  ratios <- sapply(1:length(dists), function(i) {
    median(as.numeric(as.matrix(res$dhat[[i]])/dists2[[i]]), na.rm=TRUE)
  })
  dshat <- do.call(rbind, lapply(res$cweights, diag))^2/ratios^2
  lala <- nlm(of_conv, as.numeric(rbind(res$gspace, dshat)))
  temp <- matrix(lala$estimate, n + k, p)
  xhat <- temp[1:n, ]
  dshat <- temp[-(1:n), ]
  #c(of(xhat, dshat, dists), lala$minimum)
  list(xhat = xhat, dshat = dshat, of = of, minimum = lala$minimum)
}
