

library(rgl)
library(oro.nifti)
library(prodlim)

get_fi <- function(sub, run, cope) {
  subst <- c("001","002","003","004","005","006","007","008",
             "009","010","011","012","013","014","015","016")
  paste0("model3_copes/sub", subst[sub], "/model/model003/",
         "task001_run00", run, ".feat/reg_standard/stats/cope",
         cope, ".nii.gz")
}

# get outer points from an array (done by taking dimensionwise min elts)
get_outer_points <- function(a) {
  min_inds <- apply(a, c(1, 2), function(v) min(c(Inf, which(v > 0))))
  run_inds3a <- cbind(row(min_inds)[min_inds < Inf],
                    col(min_inds)[min_inds < Inf],
                    min_inds[min_inds < Inf])
  min_inds <- apply(a, c(1, 3), function(v) min(c(Inf, which(v > 0))))
  run_inds2a <- cbind(row(min_inds)[min_inds < Inf],
                      min_inds[min_inds < Inf],
                     col(min_inds)[min_inds < Inf])
  min_inds <- apply(a, c(2, 3), function(v) min(c(Inf, which(v > 0))))
  run_inds1a <- cbind(min_inds[min_inds < Inf],
                    row(min_inds)[min_inds < Inf],
                    col(min_inds)[min_inds < Inf])
  min_inds <- apply(a, c(1, 2), function(v) max(c(-Inf, which(v > 0))))
  run_inds3b <- cbind(row(min_inds)[min_inds > -Inf],
                     col(min_inds)[min_inds > -Inf],
                     min_inds[min_inds > -Inf])
  min_inds <- apply(a, c(1, 3), function(v) max(c(-Inf, which(v > 0))))
  run_inds2b <- cbind(row(min_inds)[min_inds > -Inf],
                      min_inds[min_inds > -Inf],
                     col(min_inds)[min_inds > -Inf])
  min_inds <- apply(a, c(2, 3), function(v) max(c(-Inf, which(v > 0))))
  run_inds1b <- cbind(min_inds[min_inds > -Inf],
                     row(min_inds)[min_inds > -Inf],
                     col(min_inds)[min_inds > -Inf])
  run_inds <- rbind(run_inds1a, run_inds2a, run_inds3a,
                    run_inds1b, run_inds2b, run_inds3b)
  run_inds <- unique(run_inds)
  return(run_inds)
}

