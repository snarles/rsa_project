## get averaged brain response maps then obtain clusters

source("birds_analysis/dataloading.R")

lbda <- 20

ind <- 1
m_allz <- list()
for (ind in 1:13){
  brain <- get_brain(ind)
  phtyp <- get_typicality("ph", lbda, ind)
  pstyp <- get_typicality("ps", lbda, ind)
  m_all <- variance_brainmap(brain)
  m_allz[[ind]] <- m_all
  # m_ph <- projection_brainmap(brain, cbind(1, phtyp))
  # m_ps <- projection_brainmap(brain, cbind(1, pstyp))
  # m_both <- projection_brainmap(brain, cbind(1, pstyp, phtyp))
  # pars <- paramz[[ind]]
  # m_h <- projection_brainmap(brain, pars[, 1])
  # m_a <- projection_brainmap(brain, pars[, 2])
  # (reso <- c(phv=sum(m_ph), psv=sum(m_ps), 
  #            bv=sum(m_both), hv=sum(m_h), av=sum(m_a))/sum(m_all))

}

save(m_allz, file = "birds_analysis/m_allz.rda")


m_allz <- lapply(m_allz, log)


plot(sort(m_allz[[1]]), type = "l", main = 1)
plot(sort(m_allz[[2]]), type = "l", main = 2)


plot3d(which(m_allz[[1]] > 17, arr.ind = TRUE))
#points3d(which(m_allz[[2]] > 16, arr.ind = TRUE), col = "red")
points3d(which(m_allz[[4]] > 17, arr.ind = TRUE), col = "yellow")


plot3d(which(m_allz[[2]] > 17, arr.ind = TRUE), col = "red")
#points3d(which(m_allz[[1]] > 17, arr.ind = TRUE))
#points3d(which(m_allz[[3]] > 17, arr.ind = TRUE), col = "blue")
points3d(which(m_allz[[4]] > 17, arr.ind = TRUE), col = "yellow")


