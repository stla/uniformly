library(uniformly)
library(cgalMeshes)
library(rgl)

a <- 6; b <- 4; c <- 5

smesh <- sphereMesh(iterations = 6L)
emesh_r <- scale3d(smesh, a, b, c)
emesh <- cgalMesh$new(emesh_r)
A <- emesh$area()

capmesh <- emesh$clipToPlane(c(0, 0, -1), c(0, 0, -1), clipVolume = FALSE)
caprmesh <- capmesh$getMesh()
shade3d(caprmesh, color = "gold")
capmesh$area()

nsims <- 1e6
sims <- id(nsims, diag(c(1/a^2,1/b^2,1/c^2)), r = 1)
mean(sims[, 3] > -1) * A




nsims <- 1e6
ssims <- runif_on_sphere(nsims, 3)
allsims <- cbind(a*ssims[,1], b*ssims[,2], c*ssims[,3])
p <- sqrt((a*c*ssims[,2])^2 + (a*b*ssims[,3])^2 + (b*c*ssims[,1])^2) / b/c
e <- runif(nsims) < p 
sims <- allsims[e, ]
nsims <- nrow(sims)
summary(sims)
sum(sims[, 3] > -1) / nsims

sims <- runif_on_ellipsoid(nsims, diag(c(1/a^2,1/b^2,1/c^2)), r = 1)

nsims <- 1e7

u <- runif(nsims, -1, 1)
theta <- runif(nsims, 0, 2*pi)
x <- a * sqrt(1-u^2) * cos(theta)
y <- b * sqrt(1-u^2) * sin(theta)
z <- c * u
sims <- cbind(x,y,z)
summary(sims)
sum(sims[, 3] > -1) / nsims
