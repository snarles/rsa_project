source("birds_analysis/dataloading.R")

lambda_curvez <- list()

for (ind in 1:13){
  brain <- get_brain(ind)
  lbdas <- c(2,5,10,20,30)
  nl <- length(lbdas)
  
  mat_sub <- matrix(0, nl, 3)
  t1 <- proc.time()
  
  for (i in 1:nl) {
    lbda1 <- lbdas[i]
    lbda2 <- lbdas[i]
    phtyp <- get_typicality("ph", lbda1, ind)
    pstyp <- get_typicality("ps", lbda2, ind)
    m_ph <- projection_brainmap(brain, cbind(1, phtyp))
    m_ps <- projection_brainmap(brain, cbind(1, pstyp))
    m_both <- projection_brainmap(brain, cbind(1, pstyp, phtyp))
    mat_sub[i, ] <- c(sum(m_ph), sum(m_ps), sum(m_both))
  }
  proc.time() - t1
  
  
  lambda_curves <- cbind(lbdas, mat_sub)
  lambda_curvez[[ind]] <- lambda_curves
}