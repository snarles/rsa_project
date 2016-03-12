source("generate_copes_data.R")
source("rsa_boot_source.R")
source("ttest_source.R")
library(MASS)

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
paramz <- generate_copes_params(shared_coeffs = 0.4,
                                shared_Sigmas = 0.0,
                                mixture_param = 0.0, 
                                seed = 1, nrois = 28, 
                                nrepeats = 5
)
set.seed(321)
res <- do.call2(generate_copes_data, paramz, scale = TRUE)
zattach(res)
n <- dim(Ymat)[1]
rdat <- Ymat
## compute distance matrices
nsubs <- 16
dmats <- array(0, c(nsubs, nrois, 3))


for (sub in 1:nsubs) {
  sub.dat <- Ymat[hdat[, "sub"]==sub, ]
  sub.X <- Xmat[hdat[, "sub"]==sub, ]
  for (roi.ind in 1:nrois) {
    roi.dat <- sub.dat[, cls[, "clus"]==roi.ind]
    ##bt <- ginv(sub.X) %*% roi.dat
    ##dmat <- bt %*% t(bt)/dim(roi.dat)[2]
    dmat <- estimate_M(sub.X, roi.dat)/dim(roi.dat)[2]
    dmats[sub, roi.ind, ] <- dmat[upper.tri(dmat, diag = TRUE)]
  }
}

pvals <- array(NA, c(nrois, nrois, 3))
coord.ind <- 1; i <- 1; j <- 2
for (coord.ind in 1:3) {
  dd <- dmats[, , coord.ind]
  for (i in 1:nrois) {
    for (j in 1:nrois) {
      if (i != j) {
        tres <- t.test(dd[, i], dd[, j], paired = TRUE)
        pvals[i, j, coord.ind] <- tres$p.value
      }
    }
  }
}

image(fliplr(pvals[, , 1] < 0.05))
image(fliplr(pvals[, , 2] < 0.05))
image(fliplr(pvals[, , 3] < 0.05))

