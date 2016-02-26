####
##  Generate data with similar structure as gambling data
##  All ROIS must have same dim
####

library(pracma)

rorthog <- function(n) svd(randn(n))$u


####
#  A function with a lot of options!
#  Parameters:
#   npca : number of dimensions per roi
#   mixture_param : set to 0 to null model, 1 for alternative model
#   Wishart_mult: [1, Inf] increase to get identity (diagonal) covariances
#   Wishart_scaling: [0, Inf] increase to decrease the variance of an extra random
#    scaling factor
#  Specify:
#   'nrepeats'
#     AND
#   'params' OR 'ncopes' and 'q' (generates params matrix)
#     AND
#   'coeffs_true' and 'cls' OR  'q', nrois', 'npca',
#    'nsubjects' (generates random params) and 'mixture_param'
#     AND
#   'Sigmas_true' and 'cls' OR  'nrois', 'npca', 'nsubjects'
#     and 'Wishart_mult', and 'Wishart_scaling'
#     (uses different Sigmas for each subject)
generate_copes_params <- function(nrepeats = 3, 
                                nsubjects = NULL, ncopes = NULL, nrois = NULL,
                                npca = NULL, mixture_param = 0,
                                Wishart_mult = 2, Wishart_scaling = 100,
                                # end params for auto generation, begin custom params
                                params = NULL, cls = NULL,
                                coeffs_true = NULL, Sigmas_true = NULL,
                                ) {
  ## fill in auto-gen params using custom params
  if (!is.null(coeffs_true)) {
    nsubjects <- dim(coeffs_true)[1]
    q <- dim(coeffs_true[[1]])[1]
    p <- dim(coeffs_true[[1]])[2]
  }
  if (!is.null(Sigmas_true)) {
    nsubjects <- length(Sigmas_true)
    p <- dim(Sigmas_true)[1]
  }
  if (!is.null(params)) {
    ncopes <- dim(params)[1]
    q <- dim(params)[2]
  }
  if (!is.null(cls)) {
    nrois <- max(cls[, "clus"])
    npca <- unique(table(cls[, "clus"]))
    if (length(npca) > 1) stop("ROIs have different numbers of voxels")
    p <- dim(cls)[1]
  }
  p <- npca * nrois
  ## generate custom params if null
  if (is.null(params)) {
    params <- randn(ncopes, q)
    colnames(params) <- paste0("param", 1:q)
    rownames(params) <- paste0("cope", 1:ncopes)
  }
  if (is.null(cls)) {
    cls <- matrix(0, nrois * npca, 4)
    colnames(cls) <- c("X", "Y", "Z", "clus")
    cls[, 1:3] <- NA
    cls[,4] <- rep(1:nrois, each = npca)
    rownames(cls) <- paste0("clus", cls[,4], ".PC", rep(1:npca, nrois))
  }
  if (is.null(coeffs_true)) {
    coeffs_true <- vector(mode = "list", length = nsubjects)
    for (i in 1:nsubjects) {
      coeffs_true[[i]] <- zeros(q, p)
      common_dir <- randn(q, npca)
      for (j in 1:nrois) {
        unique_coeff <- randn(q, npca)
        cd_rot <- common_dir %*% rorthog(npca)
        coeffs_true[[i]][, cls[, "clus"]==j] <- 
          mixture_param * unique_coeff + 
          (1 - mixture_param) * cd_rot
      }
      colnames(coeffs_true[[i]]) <- rownames(cls)
      rownames(coeffs_true[[i]]) <- colnames(params)
    }
  }
  if (is.null(Sigmas_true)) {
    Sigmas_true <- vector(mode = "list", length = nsubjects)
    for (i in 1:nsubjects) {
      Sigmas_true[[i]] <- cov(randn(Wishart_mult * p)) * 
        (rgamma(1, shape = Wishart_scaling)/Wishart_scaling)
    }
  }
  list(nrepeats = nrepeats, 
       nsubjects = nsubjects, ncopes = ncopes, nrois = nrois,
       npca = npca,
       params = params, cls = cls,
       coeffs_true = coeffs_true, Sigmas_true = Sigmas_true)
}