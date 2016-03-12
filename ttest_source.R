####
##  T-test of distance matrices source code
####

## unbiased estimate of M

estimate_M <- function(X, Y) {
  B <- solve(t(X) %*% X, t(X) %*% Y)
  Yh <- X %*% B
  Xinv <- solve(t(X) %*% X)
  resid <- Y - Yh
  resD <- apply(resid, 2, var)
  B %*% t(B) - Xinv * sum(resD)
}
