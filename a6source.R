####
##  Source for analysis 6
##  Replaces role of rsa_boot_source
####

source("rsa_boot_source.R")  ## boot_sampler, sample_moments will get overwritten

## note: nX and nY args in sampler are broken
boot_sampler <- function(res) {
  p <- res$p; q <- res$q; dat <- res$dat; dat0 <- dat
  blksX <- res$blksX; blksY <- res$blksY; subsX <- res$subsX; subsY <- res$subsY
  nX <- sum(dat[, 1] == 0); nY <- sum(dat[, 1] == 1);
  nX0 <- nX; nY0 <- nY
  rawX <- dat[dat[, 1]==0, -1]
  rawY <- dat[dat[, 1]==1, -1]
  sampler <- function(nX = nX0, nY = nY0) {
    if (nX == 0 & nY == 0) {
      return(list(p = p, q = q, dat = dat, blksX = blksX, blksY = blksY, 
                  subsX = subsX, subsY = subsY))
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
    list(p = p, q = q, dat = dat, blksX = blksX, blksY = blksY, 
         subsX = subsX, subsY = subsY)
  }
  sampler
}

## computes a sample moment for each subject
sample_moments <- function(res) {
  p <- res$p; q <- res$q; dat <- res$dat
  subsX <- res$subsX; subsY <- res$subsY
  stopifnot(max(subsX) == max(subsY))
  ns <- max(subsX)
  rawX <- dat[dat[,1] == 0, -1, drop = FALSE]
  rawY <- dat[dat[,1] == 1, -1, drop = FALSE]
  Ahats <- list()
  Bhats <- list()
  M_As <- list()
  M_Bs <- list()
  for (i in 1:ns) {
    filtX <- subsX == i
    filtY <- subsY == i
    Xc <- rawX[filtX, 1:q]
    Xr <- rawX[filtX, -(1:q)]
    Yc <- rawY[filtY, 1:q]
    Yr <- rawY[filtY, -(1:q)]
    Ahat <- t(solve(t(Xc) %*% Xc, t(Xc) %*% Xr))
    Bhat <- t(solve(t(Yc) %*% Yc, t(Yc) %*% Yr))
    Xrhat <- Xc %*% t(Ahat)
    Yrhat <- Yc %*% t(Bhat)
    Xinv <- solve(t(Xc) %*% Xc)
    Yinv <- solve(t(Yc) %*% Yc)
    Xresid <- Xr - Xrhat
    Yresid <- Yr - Yrhat
    XresD <- apply(Xresid, 2, var)
    YresD <- apply(Yresid, 2, var)
    M_A <- t(Ahat) %*% Ahat - Xinv * sum(XresD)
    M_B <- t(Bhat) %*% Bhat - Yinv * sum(YresD)
    Ahats[[i]] <- Ahat; Bhats[[i]] <- Bhat; M_As[[i]] <- M_A; M_Bs[[i]] <- M_B
  }
  list(Ahats = Ahats, Bhats = Bhats, M_As = M_As, M_Bs = M_Bs)
}

## unbiased, pooled
stat.Su <- function(res) {
  mus <- sample_moments(res)
  M_As <- mus$M_As; M_Bs <- mus$M_Bs
  M_A <- Reduce(`+`, M_As)/length(M_As)
  M_B <- Reduce(`+`, M_Bs)/length(M_Bs)
  stat.S.raw <- M_A - M_B
  as.numeric(stat.S.raw)
}


jackknife <- function(res, theta) {
  theta0 <- theta(res)
  n <- dim(res$dat)[1]
  assig <- res$dat[, 1]
  blkss <- numeric(n)
  blkss[assig==0] <- res$blksX
  blkss[assig==1] <- res$blksY
  subss <- numeric(n)
  subss[assig==0] <- res$subsX
  subss[assig==1] <- res$subsY
  thetas <- list()
  for (i in 1:n) {
    res2 <- res
    res2$dat <- res$dat[-i, , drop = FALSE]
    blkss2 <- blkss[-i]; subss2 <- subss[-i]; assig2 <- assig[-i]
    res2$blksX <- blkss2[assig2==0]
    res2$blksY <- blkss2[assig2==1]
    res2$subsX <- subss2[assig2==0]
    res2$subsY <- subss2[assig2==1]
    thetas[[i]] <- theta(res2)
  }
  thetas <- do.call(cbind, thetas)
  theta_bar <- rowMeans(thetas)
  diffs <- thetas - theta_bar
  bias <- (n-1) * (theta_bar - theta0)
  sdv <- sqrt((n-1) * rowMeans((thetas - theta_bar)^2))
  skew <- rowSums(diffs^3)/(rowSums(diffs^2))^(3/2)  
  list(bias = bias, sdv = sdv, skew = skew)
}
