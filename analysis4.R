####
##  Analysis of gambling data
##  Develop bootstrap methodology
####

source("prepare_dim_reduced.R")
source("rsa_boot_source.R")

npca <- 10
p <- npca
q <- 2
res <- prepare_gambling_data(stdz_within = TRUE, npca = npca, avg_subjects = FALSE,
                             div_sqrt_p = TRUE)
zattach(res)
n <- dim(Ymat)[1]

pvals <- matrix(-1, 28, 28)
r1 <- 1; r2 <- 2
for (r1 in 1:(nrois - 1)) {
  for (r2 in (r1 + 1):nrois) {
    Y_A <- Ymat[, cls[, "clus"] == r1]
    Y_B <- Ymat[, cls[, "clus"] == r2]
    X_A <- Xmat
    X_B <- Xmat
    composite <- rbind(cbind(0, X_A, Y_A), cbind(1, X_B, Y_B))
    res <- list(p = p, q  = q, dat = composite)
    nX <- n; nY <- n
    t1 <- proc.time()
    test_res <- inverse_bca_test(res, n, n, stat.S, mc.reps = 1000)
    proc.time() - t1
    pvals[r1, r2] <- test_res
  }
}