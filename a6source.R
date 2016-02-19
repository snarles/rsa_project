####
##  Source for analysis 6
##  Replaces role of rsa_boot_source
####

source("rsa_boot_source.R")  ## boot_sampler, sample_moments will get overwritten

## note: nX and nY args in sampler are broken
boot_sampler <- function(res) {
  p <- res$p; q <- res$q; dat <- res$dat; dat0 <- dat
  blksX <- res$blksX; blksY <- res$blksY
  nX <- sum(dat[, 1] == 0); nY <- sum(dat[, 1] == 1);
  nX0 <- nX; nY0 <- nY
  rawX <- dat[dat[, 1]==0, -1]
  rawY <- dat[dat[, 1]==1, -1]
  sampler <- function(nX = nX0, nY = nY0) {
    if (nX == 0 & nY == 0) {
      return(list(p = p, q = q, dat = dat, blksX = blksX, blksY = blksY))
    }
    newX <- rawX; newY <- rawY
    for (i in 1:max(blksX)) {
      filt <- (blksX == i)
      sub <- rawX[filt, , drop = FALSE]
      newX[filt, ] <- sub[sample(dim(sub)[1], dim(sub)[1], TRUE), , drop = FALSE]
    }
    for (i in 1:max(blksY)) {
      filt <- (blksY == i)
      sub <- rawY[filt, , drop = FALSE]
      newY[filt, ] <- sub[sample(dim(sub)[1], dim(sub)[1], TRUE), , drop = FALSE]
    }
    #newX <- rawX[sample(nX0, nX, TRUE), , drop = FALSE]
    #newY <- rawY[sample(nY0, nY, TRUE), , drop = FALSE]
    dat <- rbind(cbind(0, newX), cbind(1, newY))
    list(p = p, q = q, dat = dat, blksX = blksX, blksY = blksY)
  }
  sampler
}