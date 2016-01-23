####
##  analysis of gambling data
####

library(pracma)
library(lineId)
cls <- readRDS("roi/cl_inds.rds")
params <- readRDS("roi/params.rds")
dat <- readRDS("roi/data.rds")
rdat <- dat[, -(1:3)]
hdat <- dat[, 1:3]


####
##  Classifying using models for individual subjects
####

roi_ind <- 1
sub_ind <- 7
roi_dat <- rdat[hdat[, "sub"] == sub_ind, cls[, "clus"] == roi_ind]
z_true <- hdat[hdat[, "sub"] == sub_ind, "cope"]

## strength of regularization
indiv_scale <- 0
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
res$score
