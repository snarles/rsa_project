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

scores_mat <- function(indiv_scale, quad_scale, l2p) {
  scores <- matrix(1, 16, nrois)
  
  for (roi_ind in 1:nrois) {
    for (sub_ind in 1:16) {
      roi_dat <- rdat[hdat[, "sub"] == sub_ind, cls[, "clus"] == roi_ind]
      z_true <- hdat[hdat[, "sub"] == sub_ind, "cope"]
      params2 <- cbind(scale(params), quad_scale * scale(params^2), 
                       quad_scale * scale(params[, 1] * params[, 2]),
                       indiv_scale * eye(16))
      
      ## normalize the roi_dat
      roi_dat <- scale(roi_dat, TRUE, TRUE)

      cv_scores <- numeric(3)
      # test (1, 2, 3)
      for (te_run in 1:3) {
        te_inds <- (te_run - 1)* 16 + 1:16
        X_tr <- repmat(params2, 2, 1)
        X_te <- params2
        Y_tr <- roi_dat[-te_inds, ]
        Y_te <- roi_dat[te_inds, ]
        
        res <- identification_pipeline1(X_tr, Y_tr, X_te, Y_te, 1:16,  
                                        forward_method = fit_ridge_kernel,
                                        forward_params = list(lambda = l2p),
                                        Sigma_e_method = residual_offdiag,
                                        Sigma_e_params = list(shrink = 1),
                                        backward_method = pre_mle,
                                        scoring_method = resample_misclassification,
                                        scoring_params = list(m = 2))
        cv_scores[te_run] <- res$score
      }
      scores[sub_ind, roi_ind] <- mean(cv_scores)
    }
  }
  scores
}

scores <- scores_mat(indiv_scale=1, quad_scale=0, l2p=1e6)
colMeans(scores)
mean(colMeans(scores)[1:5])
mean(colMeans(scores))

####
##  Average multiple subjects
####

avgd_scores <- function(indiv_scale, l2p) {
  av_scores <- numeric(nrois)
  for (roi_ind in 1:nrois) { 
    roi_dat <- rdat[, cls[, "clus"] == roi_ind]
    av_all <- 1/16 * repmat(eye(48), 1, 16) %*% roi_dat
    dim(roi_dat)
    dim(av_all)
    av_all <- scale(av_all, TRUE, TRUE)
    params2 <- cbind(1, scale(params),
                     indiv_scale * eye(16))
    
    cv_scores <- numeric(3)
    # test (1, 2, 3)
    for (te_run in 1:3) {
      te_inds <- (te_run - 1)* 16 + 1:16
      X_tr <- repmat(params2, 2, 1)
      X_te <- params2
      Y_tr <- av_all[-te_inds, ]
      Y_te <- av_all[te_inds, ]
      
      res <- identification_pipeline1(X_tr, Y_tr, X_te, Y_te, 1:16,  
                                      forward_method = fit_ridge_kernel,
                                      forward_params = list(lambda = l2p),
                                      Sigma_e_method = residual_offdiag,
                                      Sigma_e_params = list(shrink = 1),
                                      backward_method = pre_mle,
                                      scoring_method = resample_misclassification,
                                      scoring_params = list(m = 2))
      cv_scores[te_run] <- res$score
    }
    av_scores[roi_ind] <- mean(cv_scores)
  }
  av_scores
}

av_scores <- avgd_scores(indiv_scale = 0, l2p= 1e6)
av_scores
mean(av_scores)
(valid_rois <- which(av_scores < .5))

###
#  Conclusion: linear model does best, no need for idiosyncratic factors
###


####
##  SUFFICIENT DIMENSIONALITY REDUCTION
##  Goal: reduce the voxel dimensionality, keeping the classification performance
####

cv_score_func <- function(av_trans, l2p) {
  cv_scores <- numeric(3)
  # test (1, 2, 3)
  for (te_run in 1:3) {
    te_inds <- (te_run - 1)* 16 + 1:16
    X_tr <- repmat(params2, 2, 1)
    X_te <- params2
    Y_tr <- av_trans[-te_inds, ]
    Y_te <- av_trans[te_inds, ]
    
    res <- identification_pipeline1(X_tr, Y_tr, X_te, Y_te, 1:16,  
                                    forward_method = fit_ridge_kernel,
                                    forward_params = list(lambda = l2p),
                                    Sigma_e_method = residual_offdiag,
                                    Sigma_e_params = list(shrink = 1),
                                    backward_method = pre_mle,
                                    scoring_method = resample_misclassification,
                                    scoring_params = list(m = 2))
    cv_scores[te_run] <- res$score
  }
  mean(cv_scores)  
}

###
#  FIT THE MODEL
###

roi_ind <- 5
roi_dat <- rdat[, cls[, "clus"] == roi_ind]
av_all <- 1/16 * repmat(eye(48), 1, 16) %*% roi_dat
dim(roi_dat)
dim(av_all)
av_all <- scale(av_all, TRUE, TRUE)
params2 <- cbind(1, scale(params))
res <- svd(av_all)
av_trans <- res$u %*% diag(res$d)

cv_score_func(av_all, 1e6)
mcs <- numeric()
for (i in 10:30) {
  mcs[paste(i)] <- cv_score_func(av_trans[, 1:i], 1e6)
}
layout(matrix(1:2, 2, 1))
plot(cumsum(res$d^2), type = "l")
scatter.smooth(as.numeric(names(mcs)), mcs, xlim = c(1, 48))
fit <- loess.smooth(as.numeric(names(mcs)), mcs)
(k_chosen <- round(min.k <- fit$x[order(fit$y)[1]]))
