####
##  MAKING A DOPPELGANGER
##  Synthetic data which tries to match the properties of the real data
####

library(MASS)

seed <- 0
set.seed(0)

imp0 <- function(v) {
  v[is.na(v)] <- sample(v[!is.na(v)], sum(is.na(v)), TRUE)
  v }

dat <- readRDS("roi/data.rds")
rdat <- dat[, -(1:3)]
zpattern <- (rdat == 0)
rdat[zpattern] <- NA
rdat <- apply(rdat, 2, imp0)

mu.c <- colMeans(rdat)
mu.r <- rowMeans(rdat)
rdat.c <- t(t(rdat) - mu.c) - mu.r
Sigma0 <- cov(rdat.c)
##image(Sigma0)
##hist(diag(Sigma0))
sh <- mean(diag(Sigma0))
alpha <- 0.5
Sigma <- (1-alpha) * Sigma0 + alpha * sh

ee <- mvrnorm(n = dim(rdat)[1], mu = mu.c, Sigma = Sigma)
rdat2 <- ee + mu.r - mean(mu.r)
rdat2[zpattern] <- 0
dat2 <- cbind(dat[, 1:3], rdat2)

saveRDS(dat2, file = paste0("doppel/doppel", seed, ".rds"))
