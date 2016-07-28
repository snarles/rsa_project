library(oro.nifti)
library(AnalyzeFMRI)
library(fmri)
library(pracma)

load("birds_analysis/computed_typs.RData", verbose = TRUE)

get_brain <- function(ind) {
  fnos2 <- c("02", "03", "04", "05", "06", "07", "08", "09", "11",
             "12", "13", "14", "17")
  f <- fnos2[ind]
  fpath <- paste0("DavisPoldrack_Birds_Archive/betaseries_MNI/0", f, "_lsone_MNI.nii.gz")
  res <- oro.nifti::readNIfTI(fpath)
  res
}

get_typicality <- function(type = c("ph", "ps", "both"), lbda = 10, ind) {
  type <- type[1]
  if (type == "ph") {
    mat <- phtypz[[ind]]
    ans <- exp(-lbda * mat$rs)
  }
  if (type == "ps") {
    mat <- pstypz[[ind]]
    ans <- exp(-lbda * mat$rs)
  }
  if (type == "both") {
    mat <- phtypz[[ind]]
    ans1 <- exp(-lbda * mat$rs)
    mat <- pstypz[[ind]]
    ans2 <- exp(-lbda * mat$rs)
    ans <- cbind(ans1, ans2)
  }
  ans
}

proj_mat <- function(v) {
  if (is.null(dim(v))) v <- t(t(v))
  v %*% solve(t(v) %*% v) %*% t(v)
}

projection_brainmap <- function(brain, v, center = TRUE) {
  mat <- proj_mat(v)
  if (center) {
    ans <- apply(brain, c(1,2,3), function(w) f2(mat %*% (w- mean(w))))
  }
  else {
    ans <- apply(brain, c(1,2,3), function(w) f2(mat %*% w))
  }
  ans
}

gplot <- function(sl) {
  image(t(flipud(sl)), col = gray.colors(10))
}

f2 <- function(x, y = 0) sum((x-y)^2)
