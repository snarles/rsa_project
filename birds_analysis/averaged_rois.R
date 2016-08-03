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

##save(m_allz, file = "birds_analysis/m_allz.rda")
load("birds_analysis/m_allz.rda")

log2 <- function(v) log(pmax(v, 1))
m_allz <- lapply(m_allz, log2)


ma_filt_v <- function(v, bandw = 100) {
  trunc_len <- floor(length(v)/bandw) * bandw
  v_trunc <- v[1:trunc_len]
  vmat <- matrix(v_trunc, nrow = bandw)
  vred <- colMeans(vmat)
  vred
}

endpointz <- list()
for (indx in 1:13) {
  v0 <- sort(m_allz[[indx]])
  v <- v0[v0 > 0]
  v <- ma_filt_v(v, 50)
  sdv <- sd(v)
  xx <- 1:length(v)
  lres <- loess(v ~ xx, span = 0.05, family = "symmetric")
  #plot(v, type = "l")
  #lines(lres$fitted, col = "red")
  v <- lres$fitted
  dv <- v[-1] - v[-length(v)]
  xx <- 1:length(dv)
  lres2 <- loess(dv ~ xx, span = 0.1, family = "symmetric")
  #plot(dv, type = "l", ylim = c(-1e-2, 1e-2))
  #lines(lres2$fitted, type = "l", ylim = c(-1e-2, 1e-2), col = "red")
  dv <- lres2$fitted
  d2v <- dv[-1] - dv[-length(dv)]
  # plot(d2v, type = "l", ylim = c(-1e-4, 1e-4))
  # abline(h = 4e-6, col = "red"); abline(h = -4e-6, col = "red")
  inds <- which(abs(d2v) < 1e-5)  
  endpoints <- range(v[inds + 1])
  plot(v0, type = "l")
  abline(h = endpoints[1], col = "red")
  abline(h = endpoints[2], col = "red")
  title(indx)
  endpointz[[indx]] <- endpoints
}
thresholds <- sapply(endpointz, `[[`, 2)

# plot(sort(m_allz[[1]]), type = "l", main = 1)
# plot(sort(m_allz[[2]]), type = "l", main = 2)
# plot(sort(m_allz[[7]]), type = "l", main = 7)
# plot(sort(m_allz[[12]]), type = "l", main = 12)
# 
# 
# hist(m_allz[[1]], breaks = 30)
# plot3d(which(m_allz[[1]] > 17, arr.ind = TRUE))
# #points3d(which(m_allz[[2]] > 16, arr.ind = TRUE), col = "red")
# #points3d(which(m_allz[[4]] > 17, arr.ind = TRUE), col = "yellow")
# points3d(which(m_allz[[5]] > 18, arr.ind = TRUE), col = "yellow")
# 
# 
# plot3d(which(m_allz[[2]] > 17, arr.ind = TRUE), col = "red")
# #points3d(which(m_allz[[1]] > 17, arr.ind = TRUE))
# #points3d(which(m_allz[[3]] > 17, arr.ind = TRUE), col = "blue")
# points3d(which(m_allz[[4]] > 17, arr.ind = TRUE), col = "yellow")



indx <- 5; plot3d(which(m_allz[[indx]] > thresholds[indx], TRUE), col = "black")
indx <- 8; points3d(which(m_allz[[indx]] > thresholds[indx], TRUE), col = "red")

m_allz <- lapply(m_allz, function(v) {
  v[v < 0] <- 0
})

m_avg <- Reduce(`+`, m_allz)/13
dim(m_avg)
plot(sort(m_avg), type = "l")
plot3d(which(m_avg > 16, TRUE))
indx <- 1; points3d(which(m_allz[[indx]] > thresholds[indx], TRUE), col = "blue")
indx <- 3; points3d(which(m_allz[[indx]] > thresholds[indx], TRUE), col = "green")
indx <- 5; points3d(which(m_allz[[indx]] > thresholds[indx], TRUE), col = "red")
indx <- 7; points3d(which(m_allz[[indx]] > thresholds[indx], TRUE), col = "yellow")

plot3d(which(m_avg > 16, TRUE))
indx <- 9; points3d(which(m_allz[[indx]] > thresholds[indx], TRUE), col = "blue")
indx <- 10; points3d(which(m_allz[[indx]] > thresholds[indx], TRUE), col = "green")
indx <- 12; points3d(which(m_allz[[indx]] > thresholds[indx], TRUE), col = "red")
indx <- 13; points3d(which(m_allz[[indx]] > thresholds[indx], TRUE), col = "yellow")
