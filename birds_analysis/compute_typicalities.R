library(rgl)

fnos <- c("002", "003", "004", "005", "006", "007", "008", "009", "011",
          "012", "013", "014", "017")

fnos2 <- c("02", "03", "04", "05", "06", "07", "08", "09", "11",
           "12", "13", "14", "17")

stypz <- list()
paramz <- list()

for (f in fnos2) {
  fpath2 <- paste0("DavisPoldrack_Birds_Archive/attribute_files/", f, "_angleheight.txt")
  params <- read.table(fpath2, header = FALSE, sep = " ")
  colnames(params) <- c("angle", "height", "run")
  paramz[[f]] <- params
}

for (f in fnos) {
  fpath <- paste0("DavisPoldrack_Birds_Archive/attribute_files/", f, "_subject_typ.txt")
  styps <- read.table(fpath, header = FALSE, sep = " ")
  stypz[[f]] <- styps
}

for (ind in 1:length(fnos)) {
  pdf(paste0("birds_analysis/pplots/param", ind, ".pdf"))
  plot(paramz[[ind]][, -3])
  title(ind)
  dev.off()
}

ind <- 13
mat <- cbind(paramz[[ind]][, -3], stypz[[ind]][, 1])
plot3d(mat)


## OBSERVATIONS
##  Ind Sub  Params SubTyp  
##    1   2       A U-shaped in height
##    2   3       B flat
##    3   4       A irregular
##    4   5       C upside down flattish
##    5   6       B U-shaped in height
##    6   7       D irregular
##    7   8       C wonton: more U-shaped in angle
##    8   9       B wonton
##    9   11      C irregular: more curved in angle
##   10   12      A U-shaped in height
##   11   13      B asymmetric U-shaped in height (greater typ for short)
##   12   14      C wonton: more U-shaped in angle
##   13   17      B increasing linearly in angle, u-shaped in height


####
## COMPUTE PHYSICAL TYPICALITIES
####

allp <- do.call(rbind, paramz)
apply(allp, 2, range)

## angle from -25 to 65
## height from 57.5 to 207.5

angs <- c(-25, 65)
heits <- c(57.5, 207.5)
angRANGE <- max(angs) - min(angs)
heitRANGE <- max(heits) - min(heits)
angl_h <- mean(angs)
heit_h <- mean(heits)
umix <- c(0.25, 0.75)
lmix <- c(0.75, 0.25)

mdptz <- list()
ptypz <- list()
for (i in 1:13) {
  params <- paramz[[i]]
  nest <- 1 + (params[, 1] > angl_h) + 2 * (params[, 2] > heit_h)
  mida <- ifelse(params[, 1] > angl_h, sum(umix * angs), sum(lmix * angs))
  midh <- ifelse(params[, 2] > heit_h, sum(umix * heits), sum(lmix * heits))
  mdpts <- cbind(nest, mida, midh)
  mdptz[[i]] <- mdpts
  diFF <- params[, -3] - mdpts[, -1]
  diFF[, 1] <- diFF[, 1]/angRANGE
  diFF[, 2] <- diFF[, 2]/heitRANGE
  ptyps <- cbind(diFF, rowSums(diFF^2), exp(-10 * rowSums(diFF^2)))
  ptypz[[i]] <- ptyps
  # plot(params[, -3])
  # for (i in 1:4) {
  #   points(params[nest==i, -3], col = rainbow(4)[i])
  #   points(mida[nest==i] + 2 * rnorm(sum(nest==i)), midh[nest==i] + 2 * rnorm(sum(nest==i)), 
  #          col = rainbow(4)[i], pch = "+")
  # }
}

save(mdptz, ptypz, paramz, stypz, file = "birds_analysis/computed_typs.RData")
