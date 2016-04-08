source("source1.R")
source("visualization/flatplots.R")

cls <- readRDS(file = "roi/cl_inds.rds")

mdpts <- t(sapply(1:28, function(i) colMeans(cls[cls[, "clus"]==i, 1:3])))

bmask <- readRDS(file = 'visualization/brain_surface.rds')
clmask <- readRDS(file = 'visualization/roi_surfaces.rds')

pvalsA <- readRDS(file = "a7res.rds")

dim.ind  <- 1
pvals <- pvalsA[, , dim.ind]
nrois <- 28
rownames(pvals) <- paste0("roi", 1:nrois)
colnames(pvals) <- paste0("roi", 1:nrois)


## plot all ROIs

rcols <- rainbow(nrois)
alpha_mult <- 3
colors <- c(gray(0.5, alpha_mult * 0.3), rainbow(nrois, alpha = alpha_mult * 0.05))
rois <- c(list(bmask), clmask)

for (viewind in 1:3) {
  png(paste0('a7plots/all_rois_view', viewind, '.png'), width = 960, height = 960)
  flatplots(rois, colors, viewind = viewind)
  label_roi_plot(mdpts, viewind, trad = 10, lrad = 5, lambda = 0.002, cex = 3)
  dev.off()
}


## plot indivudal ROIS
nrois
for (roi_ind in 1:nrois) {
  pvals[roi_ind, ]
  (sig_rois <- which(pvals[roi_ind, ] < 0.05))
  ns <- length(sig_rois)
  if (ns > 0) {
    lrois <- c(roi_ind, sig_rois)
    prois <- c(list(bmask), clmask[lrois])
    colors <- c(gray(0.5, alpha_mult * 0.3), 
                rgb(c(0, rep(1, ns)), c(1, rep(0, ns)), rep(0, ns + 1), 
                    alpha = alpha_mult * 0.05))
    for (viewind in 1:3) {
      png(paste0('a7plots/roi_', roi_ind, '_view', viewind, '.png'))
      flatplots(prois, colors, viewind = viewind)
      label_roi_plot(mdpts[lrois, , drop = FALSE], labs = paste(lrois), viewind, 
                     trad = 10, lrad = 5, lambda = 0.002, cex = 3)
      title(paste('Different from roi', roi_ind, 'at alpha=0.05'))
      dev.off()
    }
  }
}
