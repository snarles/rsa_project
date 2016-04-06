

library(rgl)
library(oro.nifti)
library(prodlim)
library(igraph)

get_fi <- function(sub, run, cope) {
  subst <- c("001","002","003","004","005","006","007","008",
             "009","010","011","012","013","014","015","016")
  paste0("model3_copes/sub", subst[sub], "/model/model003/",
         "task001_run00", run, ".feat/reg_standard/stats/cope",
         cope, ".nii.gz")
}

# get outer points from an array (done by taking dimensionwise min elts)
get_outer_points <- function(a, thres = 0) {
  min_inds <- apply(a, c(1, 2), function(v) min(c(Inf, which(v > thres))))
  run_inds3a <- cbind(row(min_inds)[min_inds < Inf],
                    col(min_inds)[min_inds < Inf],
                    min_inds[min_inds < Inf])
  min_inds <- apply(a, c(1, 3), function(v) min(c(Inf, which(v > thres))))
  run_inds2a <- cbind(row(min_inds)[min_inds < Inf],
                      min_inds[min_inds < Inf],
                     col(min_inds)[min_inds < Inf])
  min_inds <- apply(a, c(2, 3), function(v) min(c(Inf, which(v > thres))))
  run_inds1a <- cbind(min_inds[min_inds < Inf],
                    row(min_inds)[min_inds < Inf],
                    col(min_inds)[min_inds < Inf])
  min_inds <- apply(a, c(1, 2), function(v) max(c(-Inf, which(v > thres))))
  run_inds3b <- cbind(row(min_inds)[min_inds > -Inf],
                     col(min_inds)[min_inds > -Inf],
                     min_inds[min_inds > -Inf])
  min_inds <- apply(a, c(1, 3), function(v) max(c(-Inf, which(v > thres))))
  run_inds2b <- cbind(row(min_inds)[min_inds > -Inf],
                      min_inds[min_inds > -Inf],
                     col(min_inds)[min_inds > -Inf])
  min_inds <- apply(a, c(2, 3), function(v) max(c(-Inf, which(v > thres))))
  run_inds1b <- cbind(min_inds[min_inds > -Inf],
                     row(min_inds)[min_inds > -Inf],
                     col(min_inds)[min_inds > -Inf])
  run_inds <- rbind(run_inds1a, run_inds2a, run_inds3a,
                    run_inds1b, run_inds2b, run_inds3b)
  run_inds <- unique(run_inds)
  return(run_inds)
}

get_cluster_inds <- function(a, thres = 1e-2) {
  raw_inds <- which(a, arr.ind = TRUE)
  dd <- as.matrix(dist(raw_inds))
  am <- (dd <= 2) + 0
  for (i in 1:10) {
    am <- am %*% am  
    am <- (am > 0) + 0    
  }
  res <- eigen(am)
  filt <- (abs(res$values) > thres)
  sum(filt)
  cls <- res$vectors[, filt]
  cls <- abs(cls) > 1e-3
  res <- apply(cls, 2, which)
  lapply(res, function(v) raw_inds[v, , drop = FALSE])
}

get_cluster_inds2 <- function(a, min.size = 5) {
  raw_inds <- which(a, arr.ind = TRUE)
  dd <- pracma::pdist(raw_inds)
  am <- (dd <= 2) + 0
  res0 <- clusters(graph.adjacency(am == 1))
  ans <- lapply(1:res0$no, function(i) which(res0$membership==i))
  ans <- ans[order(-res0$csize)]
  lans <- sapply(ans, length)
  ans <- ans[lans > min.size]
  lapply(ans, function(v) raw_inds[v, , drop = FALSE])
}
