source("source1.R")

cls <- readRDS(file = "roi/cl_inds.rds")

pvalsA <- readRDS(file = "a7res.rds")

dim.ind  <- 1
pvals <- pvalsA[, , dim.ind]
nrois <- 28
rownames(pvals) <- paste0("roi", 1:nrois)
colnames(pvals) <- paste0("roi", 1:nrois)


roi_mask <- readNIfTI('roi/zstat_mask.nii.gz')
run1_1_1 <- readNIfTI(get_fi(5, 2, 11))

## plot all ROIs
plot3d(get_outer_points(run1_1_1), aspect = FALSE, size= 1,
       xlab = "X", ylab = "Y", zlab = "Z")
for (i in 1:nrois) {
  plot3d(cls[cls[, "clus"]==i, 1:3], 
         col = rainbow(nrois)[i], size = 3, add = TRUE)
}

## plot indivudal ROIS
nrois
roi_ind <- 0

roi_ind <- roi_ind + 1
pvals[roi_ind, ]
(sig_rois <- which(pvals[roi_ind, ] < 0.05))

plot3d(get_outer_points(run1_1_1), aspect = FALSE, size= 1,
       xlab = "X", ylab = "Y", zlab = "Z")
for (i in c(roi_ind, sig_rois)) {
  cc <- grey(0.9)
  if (i == roi_ind) cc <- "green"
  if (i %in% sig_rois) cc <- "red"
  plot3d(cls[cls[, "clus"]==i, 1:3], 
         col = cc, size = 3, add = TRUE)
}
