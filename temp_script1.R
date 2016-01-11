library(oro.nifti)

## base data for clustering
roi_dat <- readNIfTI('roi/zstat.nii.gz')
image(roi_dat)

## mask
roi_mask <- readNIfTI('roi/zstat_mask.nii.gz')
image(roi_mask)

## Only uncorr results in clusters
rois <- readNIfTI('roi/thresh_zstat_uncorr.nii.gz')
table(rois)
image(rois)

