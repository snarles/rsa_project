## get averaged brain response maps then obtain clusters

source("birds_analysis/dataloading.R")

lbda <- 20

ind <- 1
for (ind in 1:13){
  brain <- get_brain(ind)
  phtyp <- get_typicality("ph", lbda, ind)
  pstyp <- get_typicality("ps", lbda, ind)
  m_all <- variance_brainmap(brain)
  m_ph <- projection_brainmap(brain, cbind(1, phtyp))
  m_ps <- projection_brainmap(brain, cbind(1, pstyp))
  m_both <- projection_brainmap(brain, cbind(1, pstyp, phtyp))
  pars <- paramz[[ind]]
  m_h <- projection_brainmap(brain, pars[, 1])
  m_a <- projection_brainmap(brain, pars[, 2])
  (reso <- c(phv=sum(m_ph), psv=sum(m_ps), 
             bv=sum(m_both), hv=sum(m_h), av=sum(m_a))/sum(m_all))

}

save(lambda_curvez, file = "birds_analysis/wbt.rda")


