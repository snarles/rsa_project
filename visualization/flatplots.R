
## rois is a list of N x 3 matrices
flatplots <- function(rois, colors, viewind = c(1, 2, 3), nsegs = 20) {
  views <- rbind(c(2, 3),c(1, 3), c(2, 1))
  viewind <- viewind[1]
  vu <- views[viewind, ]
  allrois <- do.call(rbind, rois)
  mins <- apply(allrois, 2, min)
  maxs <- apply(allrois, 2, max)
  plot(rbind(mins, maxs)[, vu], axes = FALSE, ann = FALSE, col = "white", asp = 1)
  nr <- length(rois)
  for (segi in 1:nsegs) {
    zrange <- mins[viewind] + c(segi - 1, segi) * (maxs - mins)[viewind]/nsegs
    for (roino in 1:length(rois)) {
      xx <- rois[[roino]]
      filt <- (xx[, viewind] >= zrange[1]) & (xx[, viewind] < zrange[2])
      points(xx[filt, vu, drop = FALSE], pch = ".", col = colors[roino])
    }
  }
}

## repelling labels

label_roi_plot <- function(mdpts, viewind=1, labs = paste(1:nrow(mdpts)), 
                           trad = 10, lrad = 10, lambda = 0.01, cex = 1) {
  
  
  views <- rbind(c(2, 3),c(1, 3), c(2, 1))
  viewind <- viewind[1]
  vu <- views[viewind, ]

  if (nrow(mdpts) == 1) {
    text(mdpts[, vu, drop = FALSE], labs, cex = cex)
    return()
  }
  
  
  x <- mdpts[, vu]
  x0 <- x
  of <- function(x) {
    xx <- matrix(x, ncol = 2)
    dd <- pracma::pdist(xx)
    ds <- dd[upper.tri(dd)]
    sum(pmax((trad - ds)/trad, 0)) + lambda * sum((xx - x0)^2)
  }
  
  res <- optim(as.numeric(x), of)
  ans <- matrix(res$par, ncol = 2)
  dd <- pracma::pdist(ans)
  ds <- dd[upper.tri(dd)]
  print(min(ds))  
  of(as.numeric(x))
  of(as.numeric(ans))
  ##plot(x, asp = 1)
  ds <- sqrt(rowSums((x0 - ans)^2))
  for (i in 1:nrow(x)) {
    text(ans[i, , drop = FALSE], labs[i], cex = cex)
##    text(ans[i, , drop = FALSE], labs[i], col = rcols[i])
    if (ds[i] > lrad) lines(rbind(x[i, ], ans[i, ])) 
  }
}

