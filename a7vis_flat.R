source("source1.R")
source("visualization/flatplots.R")

cls <- readRDS(file = "roi/cl_inds.rds")

mdpts <- t(sapply(1:28, function(i) colMeans(cls[cls[, "clus"]==i, 1:3])))

bmask <- readRDS(file = 'visualization/brain_surface.rds')
clmask <- readRDS(file = 'visualization/roi_surfaces.rds')

pvalsA <- readRDS(file = "a7res.rds")


pdf("a7plots/rejs.pdf")
plot(1:28, 1:28, col = "white", xlab = "roi", ylab = "roi")
for (i in 1:3) {
  text(which(pvalsA[,, i] < 0.05, TRUE), c("--", "|", "O")[i])
  text(which(pvalsA[,, i] < 0.05, TRUE), c("-", "|", "O")[i])
}
dev.off()

PV <- matrix(0, 378, 3)

agg <- (pvalsA[, , 1] < 0.05) +(pvalsA[, , 2] < 0.05) + (pvalsA[, , 3] < 0.05)
rowSums(agg, na.rm = TRUE)

for (dim.ind in 1:3) {
pvals <- pvalsA[, , dim.ind]
nrois <- 28
rownames(pvals) <- paste0("roi", 1:nrois)
colnames(pvals) <- paste0("roi", 1:nrois)

pv <- pvals[upper.tri(pvals)]
sum(pv < 0.05)
PV[, dim.ind] <- pv



## plot all ROIs

# rcols <- rainbow(nrois)
# alpha_mult <- 3
# colors <- c(gray(0.5, alpha_mult * 0.3), rainbow(nrois, alpha = alpha_mult * 0.05))
# rois <- c(list(bmask), clmask)
# 
# for (viewind in 1:3) {
#   png(paste0('a7plots/all_rois_view', viewind, '.png'), width = 960, height = 960)
#   flatplots(rois, colors, viewind = viewind)
#   label_roi_plot(mdpts, viewind, trad = 10, lrad = 5, lambda = 0.002, cex = 3)
#   dev.off()
# }


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
      png(paste0('a7plots/d_', dim.ind, 'r_', roi_ind, '_view', viewind, '.png'))
      flatplots(prois, colors, viewind = viewind)
      label_roi_plot(mdpts[lrois, , drop = FALSE], labs = paste(lrois), viewind, 
                     trad = 10, lrad = 5, lambda = 0.002, cex = 3)
      title(paste('M', c("1,1", "1,2", "1,3")[dim.ind], 'different from roi', roi_ind, 'at alpha=0.05'))
      dev.off()
    }
  }
}
}