####
##  analysis of gambling data
##  v1.0 (2016-02-01)
##  loads "data.rds", removes bad voxels, standardizes *across* subjects (within?)
##  applies PCA to each roi separately 
####

library(pracma)
library(lineId)
# 
# dfile = 'roi/data.rds'
# cfile = 'roi/cl_inds.rds'
# remove_bad = TRUE
# stdz_within = FALSE
# pca = TRUE
# npca = 10
# avg_subjects = FALSE
# stdz_params = TRUE

prepare_gambling_data <- function(dfile = 'roi/data.rds', cfile = 'roi/cl_inds.rds', remove_bad = TRUE,
                                  stdz_within = FALSE, pca = TRUE, npca = 10, avg_subjects = FALSE,
                                  stdz_params = TRUE) {
  cls <- readRDS(cfile)
  params <- readRDS("roi/params.rds")
  if (stdz_params) params <- scale(params)
  dat <- readRDS(dfile)
  if (remove_bad) {
    goodvox <- which(colSums(dat == 0) == 0)
    cls <- cls[names(goodvox)[-(1:3)], ]
    dat <- dat[, goodvox]    
  }
  rdat <- dat[, -(1:3)]
  hdat <- dat[, 1:3]
  nrois <- max(cls[, "clus"])
  sub_inds <- sort(unique(hdat[, "sub"]))
  Xmat <- repmat(params, 1, 3)
  if (stdz_within) {
    newdat <- list()
    for (sind in sub_inds) {
      subdat <- rdat[hdat[, "sub"] == sind, ]
      newdat[[paste0("S", sind)]] <- scale(subdat)
    }
    rdat <- do.call(rbind, newdat)
  } else {
    rdat <- scale(rdat)
  }
  if (avg_subjects) {
    rdat <- 1/16 * repmat(eye(48), 1, 16) %*% rdat[, ]
    Xmat <- params
  }
  if (pca) {
    newdat <- list()
    cls2 <- matrix(0, 0, 4)
    colnames(cls2) <- colnames(cls)
    if (length(npca) == 1) npca <- rep(npca, nrois)
    for (rind in 1:nrois) {
      np <- npca[rind]
      roi_dat <- rdat[, cls[, "clus"] == rind]
      res <- svd(roi_dat)
      newdat[[paste0("roi", rind)]] <- res$u %*% diag(res$d)
      temp <- cbind(NA, NA, NA, rep(rind, np))
      rownames(temp) <- paste0("clus1.PC", 1:np)
      cls2 <- rbind(cls2, temp)
    }
    rdat <- do.call(cbind, newdat)
    cls <- cls2
  }
  ans <- list(rdat = rdat, hdat = hdat, params = params, cls = cls, dat = cbind(hdat, rdat))

  ans
}