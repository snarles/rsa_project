source("birds_analysis/dataloading.R")


ind <- 3
brain <- get_brain(ind)
lbda1 <- 10
lbda2 <- 10
phtyp <- get_typicality("ph", lbda1, ind)
pstyp <- get_typicality("ps", lbda2, ind)

t1 <- proc.time()
m_0 <- projection_brainmap(brain, rep(1, length(phtyp)))
m_ph <- projection_brainmap(brain, cbind(1, phtyp))
m_ps <- projection_brainmap(brain, cbind(1, pstyp))
m_both <- projection_brainmap(brain, cbind(1, pstyp, phtyp))
proc.time() - t1

f2(m_0)


dim(m_0)

gplot(m_0[, , 40])
