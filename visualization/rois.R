cl <- readRDS("roi/cl_inds.rds")

source("visualization/source.R")

cl_pts <- list()

ncl <- max(cl[, "clus"])
i <- 1
for (i in 1:ncl) {
  xx <- cl[cl[, "clus"]==i, 1:3]
  mins <- apply(xx, 2, min) - 2
  maxs <- apply(xx, 2, max) + 2
  xx2 <- t(t(xx) - mins)
  apply(xx2, 2, max)
  dim2 <- maxs - mins
  a <- array(FALSE, dim2)
  a[xx2] <- TRUE
  # rgl::plot3d(which(a, TRUE), aspect = FALSE)
  smt <- get_surface2(a, k = 10, npoints = 1e6)
  rgl::plot3d(smt, aspect = FALSE)
  smt2 <- t(t(smt) + mins)
  # rgl::plot3d(xx, aspect = FALSE)
  # rgl::points3d(smt2, col = "red", size = 1)
  cl_pts[[i]] <- smt2
}

saveRDS(cl_pts, file = "visualization/roi_surfaces.rds")
