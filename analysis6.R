####
##  Analysis of gambling data
##  Develop bootstrap methodology
####

source("prepare_dim_reduced.R")
source("a6source.R")

npca <- 10
p <- npca
q <- 2
res <- prepare_gambling_data(dfile = "roi/data.rds",                             
                             ##res <- prepare_gambling_data(dfile = "doppel/doppel0.rds",
                             stdz_within = TRUE, npca = npca, avg_subjects = FALSE,
                             div_sqrt_p = TRUE)
zattach(res)
n <- dim(Ymat)[1]
rdat <- Ymat

## make blocks
blks <- as.numeric(as.factor(apply(hdat[, c(1, 3)], 1, function(v) paste(v, collapse = "."))))

pvals <- matrix(-1, 28, 28)
for (r1 in 1:(nrois - 1)) {
  for (r2 in (r1 + 1):nrois) {
    Y_A <- Ymat[, cls[, "clus"] == r1]
    Y_B <- Ymat[, cls[, "clus"] == r2]
    X_A <- Xmat
    X_B <- Xmat
    composite <- rbind(cbind(0, X_A, Y_A), cbind(1, X_B, Y_B))
    res <- list(p = p, q  = q, dat = composite, blksX = blks, blksY = blks, 
              subsX = hdat[, "sub"], subsY = hdat[, "sub"])
    
    # b1 <- boot_sampler(res)
    # newres <- b1()
    # sm0 <- sample_moments(res)
    # smB <- sample_moments(newres)
    (pv <- inverse_bca_test(res, n, n, stat.Su, mc.reps = 1000))
    pvals[r1, r2] <- pv
}}

saveRDS(pvals, "a6res.rds")
