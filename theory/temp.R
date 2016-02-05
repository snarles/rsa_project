####
##  Figure out how SMACOF works
####

library(smacof)
vals <- 1:10
d <- as.matrix(dist(t(t(vals))))
lala <- list(d, 4 * d)
res <- indscal(lala, ndim=1, type = "ratio")
res
res$delta
res$dhat
res$confdiss
res$conf
res$cweights

table(as.numeric(as.matrix(res$confdiss[[1]])))
