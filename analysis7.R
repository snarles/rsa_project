source("prepare_dim_reduced.R")
source("rsa_boot_source.R")
source("ttest_source.R")
library(MASS)
# res <- prepare_gambling_data(dfile = "roi/data.rds",pca = FALSE,                              
#                              stdz_within = TRUE, avg_subjects = FALSE,
#                              div_sqrt_p = FALSE)
res <- prepare_gambling_data(dfile = "roi/data.rds",pca = TRUE, npca = 10,                             
                             stdz_within = TRUE, avg_subjects = FALSE,
                             div_sqrt_p = FALSE)

zattach(res)
names(res)

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

View(pvals[, , 1])
image(fliplr(pvals[, , 1] < 0.05))
image(fliplr(pvals[, , 2] < 0.05))
image(fliplr(pvals[, , 3] < 0.05))

#saveRDS(pvals, "a7res.rds")
