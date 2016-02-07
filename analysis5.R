####
##  analysis of gambling data
##  compute separate ROI distance matrices, then apply INDSCALE
####

source("prepare_dim_reduced.R")
source("indscal_source.R")
library(pracma)
library(lineId)
library(vegan)

res <- prepare_gambling_data(dfile = "roi/data.rds",                             
                             ##res <- prepare_gambling_data(dfile = "doppel/doppel0.rds",
                             stdz_within = FALSE, pca = FALSE, avg_subjects = FALSE,
                             div_sqrt_p = FALSE)
zattach(res)
rdat <- Ymat

####
##  Functions
####
get_bscale <- function(Xmat, Ymat, lambda = 1e-3) {
  Ball <- fit_ridge_kernel(Xmat, Ymat, lambda=lambda)
  Yhat <- Xmat %*% Ball
  dim(Yhat)
  resid <- Ymat - Yhat
  resid_ses <- apply(resid, 2, sd)
  #hist(resid_ses)
  colnames(Ball) <- colnames(Ymat)
  Bscale <- t(t(Ball/resid_ses))  
}

get_bootstrap_mats <- function(replace = TRUE) {
  n <- dim(rdat)[1]
  inds <- sample(n, n, replace = replace)
  Ymat <- rdat[inds, ]
  temp <- hdat[inds, ]
  Xmat <- eye(16)[temp[, "cope"], ]
  list(Xmat = Xmat, Ymat = Ymat)
}

get_S_all <- function(Bscale)
  S_all <- as.matrix(dist(Bscale))/dim(Bscale)[2]

get_S_rois <- function(Bscale) {
  S_rois <- list()
  for (i in 1:nrois) {
    Bsub <- Bscale[, cls[, "clus"] == i]
    S_rois[[i]] <- as.matrix(dist(Bsub))/dim(Bsub)[2]  
  }
  S_rois
}

####
##  Average multiple subjects
####

zattach(get_bootstrap_mats(FALSE))
Bscale <- get_bscale(Xmat, Ymat)
S_rois <- get_S_rois(Bscale)


####
##  Bootstrap results
####

zattach(get_bootstrap_mats())
Bscale <- get_bscale(Xmat, Ymat)
S_rois <- get_S_rois(Bscale)

####
##  Apply indscal
####

res <- indscal_routine(S_rois, p=2, itmax = 4000)
plot(res$dshat)
plot(res$xhat)
resp <- vegan::procrustes(params, res$xhat)
resp
plot(resp)
resp <- vegan::procrustes(res$xhat, params)
plot(resp)

summary(lm(res$xhat[, 1] ~ params))
summary(lm(res$xhat[, 2] ~ params))

####
##  A look at the variability
####

layout(matrix(1:4, 2, 2))
for (i in 1:4) {
  zattach(get_bootstrap_mats())
  ##zattach(get_bootstrap_mats(FALSE))
  Bscale <- get_bscale(Xmat, Ymat)
  S_rois <- get_S_rois(Bscale)
  res <- indscal_routine(S_rois, p=2, verbose = TRUE, itmax = 4000)
  resp <- vegan::procrustes(params, res$xhat)
  resp
  plot(resp)
}
layout(1)

plot(procrustes(params, randn(16, 2)))
vegan::procrustes(params, randn(16, 2))
