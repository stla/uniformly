library(PlaneGeometry)
library(uniformly)

A <- rbind(c(2, 1), c(1, 1))
r <- 1
sims <- runif_on_ellipsoid(10000, A, r)
alpha1 <- 10 * pi/180
alpha2 <- 60 * pi/180

ell <- EllipticalArc$new(ell, alpha1=0, alpha2=360, degrees = TRUE)
perimeter <- ell$length()


mean(atan2(sims[,2], sims[,1]) > alpha1 & atan2(sims[,2], sims[,1]) < alpha2) * perimeter

ell <- EllipseFromCenterAndMatrix(c(0, 0), A)
arc <- EllipticalArc$new(ell, alpha1, alpha2, degrees = FALSE)
arc$length()

path <- ell$path(500L)
perim <- 0
for(i in 2L:nrow(path)){
  perim <- perim + sqrt(c(crossprod(path[i-1,]-path[i,])))
}


a <- ell$rmajor
b <- ell$rminor

x <- rnorm(100000, 0, a)
y <- rnorm(100000, 0, b)
d <- sqrt(x^2/a^2 + y^2/b^2)
sims <- cbind(x/d,y/d)

#######################

runifonellipsoid <- function(n, A, r) {
  S <- A / (r * r)
  stopifnot(isSymmetric(S))
  e <- eigen(S, symmetric = TRUE)
  if(any(e$values <= 0)) stop("`S` is not positive.")
  radiisq <- 1 / e$values
  radii <- sqrt(radiisq)
  rot <- e$vectors
  xyz <- t(vapply(radii, function(r) {
    rnorm(n, 0, r)
  }, FUN.VALUE = numeric(n))) 
  d <- sqrt(colSums(xyz*xyz / radiisq))
  sims0 <- t(xyz) / d
  t(rot %*% t(sims0))
}

sims <- runifonellipsoid(100, A, 1)

plot(NULL, xlim = c(-2, 2), ylim = c(-2, 2), asp = 1)
draw(ell)
points(sims, pch = 19)