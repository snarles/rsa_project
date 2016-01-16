

## base data for clustering
roi_dat <- readNIfTI('roi/zstat.nii.gz')
image(roi_dat)

## mask
roi_mask <- readNIfTI('roi/zstat_mask.nii.gz')
image(roi_mask)
mask_inds <- which(roi_mask == 1, arr.ind = TRUE)
mask_inds <- t(t(mask_inds) - c(27, 32, 23))
mask_inds <- mask_inds * 3
plot3d(mask_inds[mask_inds[, 1] < 0, ], aspect = FALSE)

## Only uncorr results in clusters
rois <- readNIfTI('roi/thresh_zstat_uncorr.nii.gz')
table(rois)
image(rois)

rois # 53 x 63 x 46


## data
## 16 subjects, 3 runs, 16 copes
run1_1_1 <- readNIfTI(get_fi(5, 2, 11))
run1_1_1 # 91 x 109 x 91
run_inds <- which(run1_1_1 > 0, arr.ind = TRUE)
run_inds <- t(t(run_inds) - c(46, 55, 46))
run_inds <- run_inds * 2
plot3d(run_inds[run_inds[, 1] <= 0, ], aspect = FALSE)


orthographic(run1_1_1)
orthographic(roi_dat)
image(roi_dat[, , 43], col = grey.colors(10))
image(run1_1_1[, , 85], col = grey.colors(10))

image(roi_dat[, , 3], col = grey.colors(10))
image(run1_1_1[, , 1], col = grey.colors(10))

## check symmetry

image(roi_dat[27, , ], col = grey.colors(10))

image(roi_dat[20, , ], col = grey.colors(10))
image(roi_dat[34, , ], col = grey.colors(10))

image(roi_dat[13, , ], col = grey.colors(10))
image(roi_dat[41, , ], col = grey.colors(10))

image(run1_1_1[46, , ], col = grey.colors(10))

image(run1_1_1[56, , ], col = grey.colors(10))
image(run1_1_1[36, , ], col = grey.colors(10))

image(run1_1_1[66, , ], col = grey.colors(10))
image(run1_1_1[26, , ], col = grey.colors(10))

image(run1_1_1[76, , ], col = grey.colors(10))
image(run1_1_1[16, , ], col = grey.colors(10))

image(run1_1_1[86, , ], col = grey.colors(10))
image(run1_1_1[6, , ], col = grey.colors(10))


#####
##  ALIGNMENT PLOTS
#####

mask_inds <- which(roi_mask == 1, arr.ind = TRUE)
mask_inds <- t(t(mask_inds) - c(27, 32, 23))
mask_inds <- mask_inds * 3

ex_inds <- get_outer_points(run1_1_1)
ex_inds <- t(t(ex_inds) - c(46, 55, 46))
ex_inds <- ex_inds * 2

plot3d(ex_inds[ex_inds[, 1] <= 0, ], aspect = FALSE, size= 3)
plot3d(mask_inds[mask_inds[, 1] <= 0, ], aspect = FALSE, add = TRUE, col = "green")


####
##  SEPARATE ROIS INTO CLUSTERS
####

cls <- get_cluster_inds(rois==1, 3)

####
##  APPLY ALIGNMENT TRANSFORM
####
cls2 <- lapply(cls, function(m)
  t(3/2 * (t(m) - c(27, 32, 23)) + c(46, 55, 46)))

####
##  CREATE ROI FILTER
####

rois2 <- 0 * run1_1_1
cube3 <- AlgDesign::gen.factorial(c(3, 3, 3))
for (i in 1:length(cls2)) {
  a <- cls2[[i]]
  for (j in 1:dim(a)[1]) {
    v <- floor(a[j, ])
    cv <- t(t(cube3) + v)
    rois2[cv] <- i
  }
}

plot3d(get_outer_points(run1_1_1), aspect = FALSE, size= 1,
       xlab = "X", ylab = "Y", zlab = "Z")
for (i in 1:length(cls2)) {
  plot3d(which(rois2 == i, arr.ind = TRUE), 
         col = rainbow(length(cls2))[i], size = 3, add = TRUE)
}

help(writeNIfTI)
oro.nifti::writeNIfTI(rois2, "roi/rois")

nff <- readNIfTI("roi/rois.nii.gz")
nff
table(nff)

cl_inds <- which(nff >0, arr.ind = TRUE)
cl_inds <- cbind(cl_inds, nff[nff > 0])
cl_inds <- cl_inds[order(cl_inds[, 4]), ]
View(cl_inds)
saveRDS(cl_inds, "roi/cl_inds.rds")

####
## BUILD DATA TABLE
####

ncls <- max(nff)
sub <- sample(1:16, 1)
run <- sample(1:3, 1)
cope <- sample(1:16, 1)
tabs <- list()

for (sub in 1:16) {
  for (run in 1:3) {
    for (cope in 1:16) {
      runn <- readNIfTI(get_fi(sub, run, cope))
      temp <- c(sub, run, cope, runn[cl_inds[, -4]])
      tabs <- c(tabs, list(temp))
    }
  }
}
tab <- do.call(rbind, tabs)
dim(tab)
cl_names <- apply(cl_inds[, c(4, 1:3)], 1, paste, collapse = ".")
cl_names <- paste0("clus", cl_names)
colnames(tab) <- c("sub", "run","cope", cl_names)
saveRDS(tab, "roi/data.rds")

####
## COVARIATES
####

params <- rbind(
c(13. ,   6.5),
c(13. ,  10.5),
c(13. ,  14.5),
c(13. ,  18.5),
c(21. ,   6.5),
c(21. ,  10.5),
c(21. ,  14.5),
c(21. ,  18.5),
c(29. ,   6.5),
c(29. ,  10.5),
c(29. ,  14.5),
c(29. ,  18.5),
c(37. ,   6.5),
c(37. ,  10.5),
c(37. ,  14.5),
c(37. ,  18.5))

saveRDS(params, "roi/params.rds")
