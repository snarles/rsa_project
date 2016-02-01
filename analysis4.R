####
##  Analysis of gambling data
##  Develop bootstrap methodology
####

source("prepare_dim_reduced.R")
source("rsa_boot_source.R")

npca <- 10
p <- 2
q <- npca
res <- prepare_gambling_data(stdz_within = TRUE, npca = npca, avg_subjects = FALSE)
zattach(res)

r1 <- 1; r2 <- 2
for (r1 in 1:(nrois - 1)) {
  for (r2 in (r1 + 1):nrois) {
    Y_A <- Ymat[, cls[, "clus"] == r1]
    Y_B <- Ymat[, cls[, "clus"] == r2]
    X_A <- Xmat
    X_B <- Xmat
    composite <- rbind(cbind(0, X_A, Y_A), cbind(0, X_B, Y_B))
    pack <- list(p = p, q  = q, dat = composite)
    sampler <- boot_sampler(pack)
  }
}