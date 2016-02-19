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

r1 <- 1
r2 <- 2


Y_A <- Ymat[, cls[, "clus"] == r1]
Y_B <- Ymat[, cls[, "clus"] == r2]
X_A <- Xmat
X_B <- Xmat
composite <- rbind(cbind(0, X_A, Y_A), cbind(1, X_B, Y_B))
res <- list(p = p, q  = q, dat = composite, blksX = blks, blksY = blks)
nX <- n; nY <- n
test_res <- inverse_bca_test(res, n, n, stat.Su, mc.reps = 1000)
