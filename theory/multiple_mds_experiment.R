####
##  MULTIPLE MDS EXPERIMENTS
##   minimize (over x, Sigma^k) given data D^k
##   sum_ijk |(x_i - x_j)' Sigma^k (x_i - x_j) - D_ik^k|^2
##   restrict Sigma^k to be diagonal
####



library(pracma)
library(smacof)
n <- 20
p <- 10
x <- randn(n, p)
k <- 5

ds <- abs(randn(k ,p))
ds <- t(t(ds) * c(rep(1, 2), rep(0.1, p-2)))

nms <- lapply(1:k, function(i) colSums(t(x^2) * ds[i, ]))
ips <- lapply(1:k, function(i) x %*% diag(ds[i, ]) %*% t(x))
dists <- lapply(1:k, function(i) {
  v <- t(nms[[i]])
  repmat(v, n, 1) + repmat(t(v), 1, n) - 2 * ips[[i]]
})

c(dists[[2]][5, 4], t(x[4, ] - x[5, ]) %*% diag(ds[2, ]) %*% (x[4, ] - x[5, ]))

####
##  FIND LATENT VARIABLES
####

phat <- 2
xhat <- randn(n, phat)
ds_hat <- matrix(1, k, phat)

f2 <- function(x, y = 0) sum((x - y)^2)

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

res_rnd <- nlm(of_conv, as.numeric(rbind(xhat, ds_hat)))
res_rnd$minimum

res_cht <- nlm(of_conv, as.numeric(rbind(x[, 1:2], ds[, 1:2])))
res_cht$minimum

sqrt2 <- function(v) sqrt(pmax(v, 0))
dists2 <- lapply(dists, sqrt2)

res <- indscal(dists2, ndim=2, type="ratio", verbose=TRUE)
ratios <- sapply(1:length(dists), function(i) {
  median(as.numeric(as.matrix(res$dhat[[i]])/dists2[[i]]), na.rm=TRUE)
})
names(res)
dshat <- do.call(rbind, lapply(res$cweights, diag))^2/ratios^2
of(res$gspace, dshat, dists)
lala <- nlm(of_conv, as.numeric(rbind(res$gspace, dshat)))
lala$minimum

## RESULTS

of(x, ds, dists)
of(x[, 1:2], ds[, 1:2], dists)
of(xhat, ds_hat, dists)
of(0 * xhat,0 * ds_hat, dists)
of(res$gspace, dshat, dists)

res_rnd$minimum
res_cht$minimum
lala$minimum