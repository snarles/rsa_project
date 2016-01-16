####
##  analysis of gambling data
####

cls <- readRDS("roi/cl_inds.rds")
params <- readRDS("roi/params.rds")
dat <- readRDS("roi/data.rds")
rdat <- dat[, -(1:3)]

roi_ind <- 1
roi_dat <- rdat[, cls[, "clus"] == roi_ind]

####
##  Can we cluster voxels within an roi?
####

d <- dist(t(roi_dat))
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y)

####
##  Compare subjects/copes
####

d <- dist(roi_dat)
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

## subjects

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",	type="n")
cols <- rainbow(16)
for (i in 1:16) {
  points(cbind(x, y)[dat[, "sub"]==i, ], col = cols[i])
}

## copes

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",  type="n")
cols <- rainbow(16)
for (i in 1:16) {
  points(cbind(x, y)[dat[, "cope"]==i, ], col = cols[i])
}

####
##  ANOVA analysis
####

factors <- data.frame(sub = as.factor(dat[, 1]),
                      run = as.factor(dat[, 2]),
                      cope = as.factor(dat[, 3]))

anovas <- array(0, dim = c(33, 4, dim(roi_dat)[2]))
for (i in 1:dim(roi_dat)[2]) {
  temp <- lm(y ~ sub + run + cope, 
             data = data.frame(y = roi_dat[, i], factors))
  mm <- coef(summary(temp))
  anovas[, , i] <- mm
}

stuff <- (t(anovas[, 1, ]))
colnames(stuff) <- rownames(mm)
View(stuff)
matplot(t(anovas[2:16, 1, ]), type = "l")
matplot(t(anovas[17:18, 1, ]), type = "l")
matplot(t(anovas[19:33, 1, ]), type = "l")

