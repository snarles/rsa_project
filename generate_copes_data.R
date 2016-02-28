####
##  Generate data with similar structure as gambling data
##  All ROIS must have same dim
####

library(pracma)
library(MASS)

rorthog <- function(n) svd(randn(n))$u


####
#  Generate random parameters
#  Arguments:
#   npca : number of dimensions per roi
#   mixture_param : set to 0 to null model, 1 for alternative model
#   Wishart_mult: [1, Inf] increase to get identity (diagonal) covariances
#   Wishart_scaling: [0, Inf] increase to decrease the variance of an extra random
#    scaling factor
#   shared_coefs: [0, 1], 1=all subjects have the same coefficients
#   shared_Sigmas: [0, 1], 1 = all subjects have the same Sigma matrix
#   Sigma_dependence: [0, 1], 0 = all ROIs have indepedent noise
generate_copes_params <- function(nrepeats = 3, 
                                nsubjects = 16, ncopes = 16, dimension = 2, nrois = 5,
                                npca = 3, mixture_param = 0,
                                Wishart_mult = 2, Wishart_scaling = 100,
                                shared_coeffs = 1, shared_Sigmas = 1,
                                Sigma_dependence = 0, seed = NULL, ...
                                ) {
  if (!is.null(seed)) set.seed(seed)
  p <- nrois * npca
  params <- randn(ncopes, dimension)
  colnames(params) <- paste0("param", 1:dimension)
  rownames(params) <- paste0("cope", 1:ncopes)
  cls <- matrix(0, nrois * npca, 4)
  colnames(cls) <- c("X", "Y", "Z", "clus")
  cls[, 1:3] <- NA
  cls[,4] <- rep(1:nrois, each = npca)
  rownames(cls) <- paste0("clus", cls[,4], ".PC", rep(1:npca, nrois))
  ## construct true coefficients
  coeffs_true <- vector(mode = "list", length = nsubjects + 1)
  for (i in 1:(nsubjects + 1)) {
    coeffs_true[[i]] <- zeros(dimension, p)
    common_dir <- randn(dimension, npca)
    for (j in 1:nrois) {
      unique_coeff <- randn(dimension, npca)
      cd_rot <- common_dir %*% rorthog(npca)
      coeffs_true[[i]][, cls[, "clus"]==j] <- 
        mixture_param * unique_coeff + 
        (1 - mixture_param) * cd_rot
    }
    colnames(coeffs_true[[i]]) <- rownames(cls)
    rownames(coeffs_true[[i]]) <- colnames(params)
  }
  for (i in 2:(nsubjects + 1)) 
    coeffs_true[[i]] <- shared_coeffs * coeffs_true[[1]] + (1-shared_coeffs) * coeffs_true[[i]]
  coeffs_true <- coeffs_true[-1]
  ## construct the Sigmas
  roi_mask <- (1 - Sigma_dependence) * eye(nrois) %x% ones(npca) + Sigma_dependence
  Sigmas_true <- vector(mode = "list", length = nsubjects + 1)
  for (i in 1:(nsubjects + 1)) {
    Sigmas_true[[i]] <- cov(randn(Wishart_mult * p, p)) * roi_mask * 
      (rgamma(1, shape = Wishart_scaling)/Wishart_scaling)
    colnames(Sigmas_true[[i]]) <- rownames(cls)
    rownames(Sigmas_true[[i]]) <- rownames(cls)
  }
  for (i in 2:(nsubjects + 1)) 
    Sigmas_true[[i]] <- shared_Sigmas * Sigmas_true[[1]] + (1-shared_Sigmas) * Sigmas_true[[i]]
  Sigmas_true <- Sigmas_true[-1]
  list(nrepeats = nrepeats, dimension = dimension,
       nsubjects = nsubjects, ncopes = ncopes, nrois = nrois,
       npca = npca,
       params = params, cls = cls,
       coeffs_true = coeffs_true, Sigmas_true = Sigmas_true)
}

generate_copes_data <- function(nrepeats, nsubjects, ncopes, params, coeffs_true, cls,
                                Sigmas_true, ...) {
  ansY <- list()
  p <- dim(coeffs_true[[1]])[2]
  Xindiv <- repmat(params, nrepeats, 1)
  hdat <- cbind(sub = rep(1:nsubjects, each = ncopes * nrepeats),
                run = rep(rep(1:nrepeats, each = ncopes), nsubjects),
                cope = rep(1:ncopes, nrepeats * nsubjects))
  for (ii in 1:nsubjects) {
    noise <- mvrnorm(ncopes * nrepeats, rep(0, p), Sigmas_true[[ii]])
    mu <- Xindiv %*% coeffs_true[[ii]]
    ansY[[ii]] <- mu + noise
  }
  Ymat <- do.call(rbind, ansY)
  Xmat <- repmat(Xindiv, nsubjects, 1)
  list(Ymat = Ymat, Xmat = Xmat, hdat = hdat, cls = cls, nrois = nrois,
       params = params, dat = cbind(hdat, Ymat), 
       coeffs_true = coeffs_true, Sigmas_true = Sigmas_true)
}

lineId::zattach(formals(generate_copes_params))
paramz <- generate_copes_params()
lala <- do.call(generate_copes_data, paramz)
names(lala)
zattach(lala)
