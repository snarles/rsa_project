library(rgl)
n <- 3 ## try n = 20 and plot it in a large window, it'll blow your mind
ps <- which(array(TRUE, dim = rep(n * 2 + 1, 3)), arr.ind=TRUE)
ps <- (ps - n - 1)/n
ps <- ps[rowSums(ps^2) <=1, ]
plot3d(ps, size=20/n, axes = FALSE, 
       ann = FALSE, xlab = "", ylab = "", zlab = "")
