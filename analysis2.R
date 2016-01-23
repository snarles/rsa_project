####
##  analysis of gambling data
####

library(pracma)
library(lineId)
cls <- readRDS("roi/cl_inds.rds")
params <- readRDS("roi/params.rds")
dat <- readRDS("roi/data.rds")

## delete bad voxels
goodvox <- which(colSums(dat == 0) == 0)
cls <- cls[names(goodvox)[-(1:3)], ]
dat <- dat[, goodvox]

rdat <- dat[, -(1:3)]
hdat <- dat[, 1:3]
nrois <- max(cls[, "clus"])



####
##  Classifying using models for individual subjects
####

residual_offdiag2 <- function (X, Y, B, filt = rep(TRUE, dim(Y)[2]), shrink = 0.5, 
          mc.cores = 1, ...) 
{
  Y <- Y[, filt]
  Yh <- paramultiply(X, B, mc.cores)
  resids <- Y - Yh
  (1 - shrink) * cov(resids) + shrink * diag(diag(cov(resids)) + 1e-5)
}

scores <- matrix(1, 16, nrois)

for (roi_ind in 1:nrois) {
  for (sub_ind in 1:16) {
    roi_dat <- rdat[hdat[, "sub"] == sub_ind, cls[, "clus"] == roi_ind]
    z_true <- hdat[hdat[, "sub"] == sub_ind, "cope"]
    
    ## strength of regularization
    indiv_scale <- 10
    params2 <- cbind(scale(params), indiv_scale * eye(16))
    
    ## normalize the roi_dat
    roi_dat <- scale(roi_dat, TRUE, TRUE)
    
    # test (1, 2, 3)
    te_run <- 1
    te_inds <- (te_run - 1)* 16 + 1:16
    X_tr <- repmat(params2, 2, 1)
    X_te <- params2
    Y_tr <- roi_dat[-te_inds, ]
    Y_te <- roi_dat[te_inds, ]
    
    res <- identification_pipeline1(X_tr, Y_tr, X_te, Y_te, 1:16,  
                                    forward_method = fit_ridge_kernel,
                                    forward_params = list(lambda = 1e-3),
                                    Sigma_e_method = residual_offdiag,
                                    Sigma_e_params = list(shrink = 1),
                                    backward_method = pre_mle,
                                    scoring_method = resample_misclassification,
                                    scoring_params = list(m = 2))
    scores[sub_ind, roi_ind] <- res$score
  }
}
