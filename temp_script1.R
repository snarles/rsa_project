

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

rois # 53 x 63 x 46


## data
## 16 subjects, 3 runs, 16 copes
run1_1_1 <- readNIfTI(get_fi(5, 2, 11))
run1_1_1 # 91 x 109 x 91

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
