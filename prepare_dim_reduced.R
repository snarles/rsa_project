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
# div_sqrt_p = TRUE

####
##  TODO
##   * Make the code fully compatible with subsampling schemes
##   * Implement methods for choosing # of PCA components
####

#' Arguments:
#' @param dfile The file containing voxel data
#' @param cfile File containing voxel XYZ indices and cluster index
#' @param remove_bad there are some voxels with 0.00 (equivalent of NA?): remove them
#' @param stdz_within Standardize within each subject?
#' @param pca Replace voxels with principal coordinates?
#' @param npca Number of pca components for each roi, single number if same for all rois
#' @param avg_subjects Replace all subjects with one "averaged" subject
#' @param stz_params Standardize the X matrix
#' @param div_sqrt_p Normalize by number of voxels (so we can compare average distance matrices)

prepare_gambling_data <- function(dfile = 'roi/data.rds', cfile = 'roi/cl_inds.rds', remove_bad = TRUE,
                                  stdz_within = FALSE, pca = TRUE, npca = 10, avg_subjects = FALSE,
                                  stdz_params = TRUE, div_sqrt_p = TRUE) {
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
  Xmat <- params[hdat[, "cope"], ]
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
    ## INCOMPLETE -- this might not work crrectly if subjects are subsampled
    rdat <- 1/16 * repmat(eye(48), 1, 16) %*% rdat[, ]
    Xmat <- params
  }
  if (div_sqrt_p) {
    newdat <- list()
    for (rind in 1:nrois) {
      roi_dat <- rdat[, cls[, "clus"] == rind]
      newdat[[paste0("roi", rind)]] <- rdat/sqrt(dim(rdat)[2])
    }
    rdat <- do.call(cbind, newdat)
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
  ans <- list(Ymat = rdat, hdat = hdat, params = params, cls = cls, dat = cbind(hdat, rdat),
              Xmat = Xmat, nrois = nrois)

  ans
}