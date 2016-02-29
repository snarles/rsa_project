####
##  Analysis of gambling data
##  Develop bootstrap methodology
####

#source("prepare_dim_reduced.R")
source("a6source.R")
source("generate_copes_data.R")
library(lineId)


paramz <- generate_copes_params(shared_coeffs = 0.4,
                                shared_Sigmas = 0.8,
                                mixture_param = 0.1, 
                                seed = 1,
                                nrepeats = 20
                                )
set.seed(321)
res <- do.call(generate_copes_data, paramz)
zattach(res)
n <- dim(Ymat)[1]
rdat <- Ymat

## make blocks
blks <- as.numeric(as.factor(apply(hdat[, c(1, 3)], 1, function(v) paste(v, collapse = "."))))
r1 <- 1; r2 <- 2
# pvals <- matrix(-1, 28, 28)
# for (r1 in 1:(nrois - 1)) {
#   for (r2 in (r1 + 1):nrois) {
    Y_A <- Ymat[, cls[, "clus"] == r1]
    Y_B <- Ymat[, cls[, "clus"] == r2]
    X_A <- Xmat
    X_B <- Xmat
    composite <- rbind(cbind(0, X_A, Y_A), cbind(1, X_B, Y_B))
    res <- list(p = dim(Y_A)[2], q  = dim(X_A)[2], dat = composite, blksX = blks, blksY = blks, 
              subsX = hdat[, "sub"], subsY = hdat[, "sub"])
    
    b1 <- boot_sampler(res)
    newres <- b1()
    sm0 <- sample_moments(res)
    smB <- sample_moments(newres)
    stat.Su(res)
    stat.Su(newres)
    #pvals[r1, r2] <- 
      inverse_bca_test(res, n, n, stat.Su, mc.reps = 1000)
# }}

#saveRDS(pvals, "a6res.rds")
