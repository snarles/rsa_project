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


infl_point_lvs <- function(v) {
  v <- sort(m_allz[[1]])
  v <- v[v > -Inf]
  plot(v, type = "l")
  dv <- v[-1] - v[-length(v)]
  plot(dv, type = "l", ylim = c(-1e-2, 1e-2))
  d2v <- dv[-1] - dv[-length(dv)]
  plot(d2v, type = "l", ylim = c(-1e-3, 1e-3))
  abline(h = 1e-4, col = "red")
  
}

plot(sort(m_allz[[1]]), type = "l", main = 1)
plot(sort(m_allz[[2]]), type = "l", main = 2)
plot(sort(m_allz[[7]]), type = "l", main = 7)
plot(sort(m_allz[[12]]), type = "l", main = 12)


hist(m_allz[[1]], breaks = 30)
plot3d(which(m_allz[[1]] > 17, arr.ind = TRUE))
#points3d(which(m_allz[[2]] > 16, arr.ind = TRUE), col = "red")
#points3d(which(m_allz[[4]] > 17, arr.ind = TRUE), col = "yellow")
points3d(which(m_allz[[5]] > 18, arr.ind = TRUE), col = "yellow")


plot3d(which(m_allz[[2]] > 17, arr.ind = TRUE), col = "red")
#points3d(which(m_allz[[1]] > 17, arr.ind = TRUE))
#points3d(which(m_allz[[3]] > 17, arr.ind = TRUE), col = "blue")
points3d(which(m_allz[[4]] > 17, arr.ind = TRUE), col = "yellow")


