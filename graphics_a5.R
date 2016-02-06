####
##  Graphics for analysis 5
####

temp <- apply(params, 2, function(v) as.numeric(as.factor(v)))
dm <- as.matrix(dist(temp, "manhattan"))
adjmat <- (dm == 1)
image(adjmat)
pairz <- which(adjmat, TRUE)
dist0 <- as.matrix(dist(params))

delta <- S_rois[[1]]
delta <- dist0

## initialization
res <- mds(delta, ndim=15)
xhat <- res$conf[, 1:2]
delta <- as.matrix(dist(res$conf))
resp <- procrustes(xhat, params)
resp
plot(resp)
names(resp)
resp$Yrot
resp$rotation
resp$translation
resp$scale
xhat
newparams <- t(
  t(params %*% resp$rotation * resp$scale) + 
    as.numeric(resp$translation))
plot(newparams)
points(xhat, col = "red")

####
## distance-matching flow
####

plot(as.matrix(dist(newparams)), as.matrix(dist(xhat)))
plot(as.matrix(delta), as.matrix(dist(res$conf)))
plot(as.matrix(dist(newparams)), delta)

plot(params)

dist_flow <- function(x0, delta, eps = 1e-2, nits = 100, plot = FALSE) {
  if (plot) {
    plot(3 * x0, col = "white"); points(x0)
  }
  n <- dim(delta)[1]
  p <- dim(x0)[2]  
  of <- function(x) {
    conf <- matrix(x, n, p)
    diff <- delta - as.matrix(dist(conf))
    sum(diff^2)
  }
  for (i in 1:nits) {
    x0 <- x0 - eps * numDeriv::grad(of, x0)
    if (plot) {
      points(x0, col = rainbow(nits)[i], pch = ".")
    }
  }
  x0
}

x0 <- newparams
x0 <- dist_flow(x0, delta, eps = 1e-3, nits = 10, plot = TRUE)

x1 <- x0 %*% t(resp$rotation)
plot(x1)
for (pp in apply(pairz, 1, list)) {
  lines(x1[pp[[1]], ])
}

####
##  Do this for all rois
####

nits <- 20
for (roi_ind in 1:28) {
  delta <- S_rois[[roi_ind]]
  res <- mds(delta, ndim=15)
  xhat <- res$conf[, 1:2]
  delta <- as.matrix(dist(res$conf))
  resp <- procrustes(xhat, params)
  newparams <- t(
    t(params %*% resp$rotation * resp$scale) + 
      as.numeric(resp$translation))
  png(paste0("temp_plots/roi_dist_", roi_ind, ".png"))  
  layout(matrix(c(1,3, 2, 4), 2, 2))
  for (i in 1:4) {
    x0 <- newparams    
    x0 <- dist_flow(x0, delta, eps = 1e-3, nits = i * nits, plot = FALSE)
    x1 <- x0 %*% t(resp$rotation)
    plot(x1)
    for (pp in apply(pairz, 1, list)) {
      lines(x1[pp[[1]], ])
    }
    title(paste("ROI", roi_ind, "nits = ", i * nits))
  }
  dev.off()
}

