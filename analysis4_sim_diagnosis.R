####
##  Analysis of gambling data
##  Investigate phenomena
####

source("prepare_dim_reduced.R")
source("rsa_boot_source.R")

npca <- 10
p <- npca
q <- 2
##res <- prepare_gambling_data(dfile = "roi/data.rds",                             
res_div <- prepare_gambling_data(dfile = "doppel/doppel0.rds",
                             stdz_within = FALSE, npca = npca, avg_subjects = FALSE,
                             div_sqrt_p = TRUE)
res_no <- prepare_gambling_data(dfile = "doppel/doppel0.rds",
                                 stdz_within = FALSE, npca = npca, avg_subjects = FALSE,
                                 div_sqrt_p = FALSE)

zattach(res_div)
zattach(res_no)

n <- dim(Ymat)[1]
r1 <- 1; r2 <- 2
Y_A <- Ymat[, cls[, "clus"] == r1]
Y_B <- Ymat[, cls[, "clus"] == r2]
X_A <- Xmat
X_B <- Xmat
composite <- rbind(cbind(0, X_A, Y_A), cbind(1, X_B, Y_B))
res <- list(p = p, q  = q, dat = composite)
(sm <- stat.S(res))

theta <- stat.S
sps <- sampling_dist(boot_sampler(res), stat.S, n, n, 100, TRUE)
hist(sps[1, ])
hist(sps[2, ])
hist(sps[3, ])
hist(sps[4, ])

mc.reps <- 1000
jj <- jackknife(res, theta)
acc <- jj$skew/6
boot_s1 <- boot_sampler(res)
thetas <- sampling_dist(boot_s1, theta, nX, nY, mc.reps, TRUE)
theta0 <- theta(res)
z0s <- sapply(1:length(theta0), function(i) {
  qnorm(sum(thetas[i, ] < theta0[i])/mc.reps)
})
pvs <- sapply(1:length(theta0), function(i) {
  v <- thetas[i, ]
  p1 <- (sum(v > 0) + 1)/(mc.reps + 1)
  p2 <- (sum(v < 0) + 1)/(mc.reps + 1)
  p1_mod <- inv_alpha_bca(p1, z0s[i], acc[i])
  p2_mod <- inv_alpha_bca(p2, -z0s[i], acc[i])
  2 * pmin(p1_mod, p2_mod)
})
min(pvs) * length(pvs)


nX <- n; nY <- n
t1 <- proc.time()
(test_res <- inverse_bca_test(res, n, n, stat.S, mc.reps = 1000))
proc.time() - t1
