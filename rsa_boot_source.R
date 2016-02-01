####
##  Bootstrapping testing of RSA// procrustes regression
####

library(pracma)
library(MASS)

f2 <- function(x, y = 0) sum((x-y)^2)

## DATA MODEL

regression_data_model_ <- function(A, B, SigmaX, SigmaY) {
  p <- dim(A)[1]; q <- dim(A)[2]
  sampler <- function(nX, nY) {
    noise_mult <- 1
    if (nX == 0 && nY == 0) {
      nX <- 1; nY <- 1
      noise_mult <- 0
    }
    rawXc <- mvrnorm(nX, rep(0, q), SigmaX)
    rawXr <- rawXc %*% t(A) + randn(nX, p)
    rawYc <- mvrnorm(nY, rep(0, q), SigmaY)
    rawYr <- rawYc %*% t(B) + randn(nY, p)    
    dat <- rbind(cbind(0, rawXc, rawXr), cbind(1, rawYc, rawYr))
    list(p = p, q  = q, dat = dat)
  }
  # force eval
  lala <- sampler(10, 10)
  sampler
}

sample_moments <- function(res) {
  p <- res$p; q <- res$q; dat <- res$dat
  rawX <- dat[dat[,1] == 0, -1, drop = FALSE]
  rawY <- dat[dat[,1] == 1, -1, drop = FALSE]
  Xc <- rawX[, 1:q]
  Xr <- rawX[, -(1:q)]
  Yc <- rawY[, 1:q]
  Yr <- rawY[, -(1:q)]
  Ahat <- t(solve(t(Xc) %*% Xc, t(Xc) %*% Xr))
  Bhat <- t(solve(t(Yc) %*% Yc, t(Yc) %*% Yr))
  list(Ahat = Ahat, Bhat = Bhat)
}

## TEST STATISTICS

stat.T <- function(res) {
  mus <- sample_moments(res)
  res <- svd(mus$Ahat %*% t(mus$Bhat))
  R <- res$v %*% t(res$u)
  stat.T.raw <- R %*% mus$Ahat - mus$Bhat
  as.numeric(stat.T.raw)
}

stat.S <- function(res) {
  mus <- sample_moments(res)
  stat.S.raw <- t(mus$Ahat) %*% mus$Ahat - t(mus$Bhat) %*% mus$Bhat
  as.numeric(stat.S.raw)
}

nm_ <- function(theta) {
  newf <- function(res) {
    f2(theta(res))
  }
  newf
}

stat.T.nm <- nm_(stat.T)
stat.S.nm <- nm_(stat.S)

## BOOTSTRAP

sampling_dist <- function(sampler, theta, nX, nY, mc.reps = 1e3, samples = FALSE) {
  # true params
  dat0 <- sampler(0, 0)
  theta0 <- theta(dat0)
  # get population
  thetas <- lapply(1:mc.reps, function(i) theta(sampler(nX, nY)))    
  thetas <- do.call(cbind, thetas)
  if (samples) {
    return(thetas)
  }
  diffs <- thetas - theta0
  bias <- rowMeans(diffs)
  sdv <- apply(diffs, 1, sd)
  skew <- rowSums(diffs^3)/(rowSums(diffs^2))^(3/2)
  list(bias = bias, sdv = sdv, skew = skew)
}

jackknife <- function(res, theta) {
  theta0 <- theta(res)
  n <- dim(res$dat)[1]
  thetas <- list()
  for (i in 1:n) {
    dat2 <- res$dat[-i, , drop = FALSE]
    thetas[[i]] <- theta(list(p = res$p, q = res$q, dat = dat2))
  }
  thetas <- do.call(cbind, thetas)
  theta_bar <- rowMeans(thetas)
  diffs <- thetas - theta_bar
  bias <- (n-1) * (theta_bar - theta0)
  sdv <- sqrt((n-1) * rowMeans((thetas - theta_bar)^2))
  skew <- rowSums(diffs^3)/(rowSums(diffs^2))^(3/2)  
  list(bias = bias, sdv = sdv, skew = skew)
}

boot_sampler <- function(res) {
  p <- res$p; q <- res$q; dat <- res$dat; dat0 <- dat
  nX <- sum(dat[, 1] == 0); nY <- sum(dat[, 1] == 1);
  nX0 <- nX; nY0 <- nY
  rawX <- dat[dat[, 1]==0, -1]
  rawY <- dat[dat[, 1]==1, -1]
  mus <- sample_moments(res)
  ress <- svd(mus$Ahat %*% t(mus$Bhat))
  R <- ress$v %*% t(ress$u)
  sampler <- function(nX = nX0, nY = nY0) {
    if (nX == 0 & nY == 0) {
      return(list(p = p, q = q, dat = dat0))
    }
    newX <- rawX[sample(nX0, nX, TRUE), , drop = FALSE]
    newY <- rawY[sample(nY0, nY, TRUE), , drop = FALSE]
    dat <- rbind(cbind(0, newX), cbind(1, newY))
    list(p = p, q = q, dat = dat)
  }
  sampler
}

inv_alpha_bca <- function(alpha_bca, z0, acc) {
  if (alpha_bca == 1) {
    return(1)
  }
  za <- qnorm(alpha_bca)
  pnorm(
    ((1 - acc * z0) * (za - z0) - z0)/
      (1 + acc * (za - z0))
  )
}

inverse_bca_test <- function(res, nX, nY, theta, mc.reps = 1000) {
  jj <- jackknife(res, theta)
  acc <- jj$skew/6
  boot_s1 <- boot_sampler(res)
  thetas <- sampling_dist(boot_s1, theta, nX, nY, mc.reps, TRUE)
  theta0 <- theta(res)
  z0s <- sapply(1:length(theta0), function(i) {
    qnorm(sum(thetas[i, ] < theta0[i])/mc.reps)
  })
  pvs <- sapply(1:length(theta0), function(i) {
    v <- thetas[i, ]
    p1 <- (sum(v > 0) + 1)/(mc.reps + 1)
    p2 <- (sum(v < 0) + 1)/(mc.reps + 1)
    p1_mod <- inv_alpha_bca(p1, z0s[i], acc[i])
    p2_mod <- inv_alpha_bca(p2, -z0s[i], acc[i])
    2 * pmin(p1_mod, p2_mod)
  })
  min(pvs) * length(pvs)
}

####
##  Demo
####

# p <- 5
# q <- 2
# SigmaX <- 3 * cov(randn(2*q, q)); SigmaY <- 3 * cov(randn(2 * q, q))
# A_0 <- randn(p, q); B_0 <- svd(randn(p, p))$u %*% A_0
# h0_small <- regression_data_model_(A_0, B_0, SigmaX, SigmaY)
# 
# B_1 <- randn(p, q)
# h1_small <- regression_data_model_(A_0, B_1, SigmaX, SigmaY)
# 
# dat <- h0_small(200, 200)
# mus <- sample_moments(dat)
# c(f2(mus$Ahat, A_0), f2(mus$Bhat, B_0))
# 
# 
# nX <- 100; nY <- 100; mc.reps = 1000
# res0 <- h0_small(nX, nY)
# res1 <- h1_small(nX, nY)
# 
# c(inverse_bca_test(res0, stat.T, mc.reps), inverse_bca_test(res1, stat.T, mc.reps))
# c(inverse_bca_test(res0, stat.T.nm, mc.reps), inverse_bca_test(res1, stat.T.nm, mc.reps))
# c(inverse_bca_test(res0, stat.S, mc.reps), inverse_bca_test(res1, stat.S, mc.reps))
# c(inverse_bca_test(res0, stat.S.nm, mc.reps), inverse_bca_test(res1, stat.S.nm, mc.reps))
