## get averaged brain response maps then obtain clusters

source("birds_analysis/dataloading.R")

lbda <- 20

ind <- 1
## ind <- 2
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


##########

hist(m_all)
max(m_all)
## [1] 807531741
max(log(m_all))
## 20.50949
hist(log(m_all))
which(m_all == max(m_all), arr.ind = TRUE)
## dim1 dim2 dim3
## [1,]   38   49   12
plot(brain[38, 49, 12, ])
plot(brain[38, 49, 11, ])
plot(brain[38, 49, 10, ])
plot(brain[38, 49, 12, ])
plot(brain[38, 44, 12, ])
plot(sort(log(m_all)))
gplot((log(m_all) > 16)[, , 20])
gplot((log(m_all) > 16)[, , 30])
gplot((log(m_all) > 16)[, , 40])


gplot((log(m_all) > 15.2)[, , 20])
gplot((log(m_all) > 15.2)[, , 30])
gplot((log(m_all) > 15.2)[, , 40])

library(rgl)
plot3d(which((log(m_all) > 15.9), arr.ind = TRUE))

plot(sort(log(m_both)))
gplot((log(m_both) > 12.5)[, , 20])
gplot((log(m_both) > 11)[, , 30])
gplot((log(m_both) > 11)[, , 40])

plot3d(which((log(m_all) > 15.9), arr.ind = TRUE))
points3d(which((log(m_both) > 12.5), arr.ind = TRUE), add = TRUE, col = "red")
