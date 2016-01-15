


library(oro.nifti)


get_fi <- function(sub, run, cope) {
  subst <- c("001","002","003","004","005","006","007","008",
             "009","010","011","012","013","014","015","016")
  paste0("model3_copes/sub", subst[sub], "/model/model003/",
         "task001_run00", run, ".feat/reg_standard/stats/cope",
         cope, ".nii.gz")
}
