library(rgl)
ps <- which(array(TRUE, dim = rep(41, 3)), arr.ind=TRUE)
ps <- (ps - 21)/20
ps <- ps[rowSums(ps^2) <=1, ]
plot3d(ps, size=1, axes = FALSE, 
       ann = FALSE, xlab = "", ylab = "", zlab = "")
