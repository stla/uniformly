
A <- rbind(c(2, 1), c(1, 1))
r <- 2



library(uniformly)
set.seed(666L)
# ellipse parameters
A <- rbind(c(2, 1), c(1, 1))
r <- 2
# plot the ellipse
x1 <- seq(-2.5, 2.5, length.out = 100)
x2 <- seq(-3, 3, length.out = 100)
z <- outer(
  x1, x2, FUN = Vectorize(function(x1, x2) t(c(x1, x2)) %*% A %*% c(x1, x2))
)
contour(x1, x2, z, nlevels = 1, levels = r^2, asp = 1, drawlabels = FALSE)
# simulations on the perimeter
sims <- runif_on_ellipsoid(60, A, r)
points(sims, pch = 19, col = "blue")
# simulations in the area
sims <- runif_in_ellipsoid(100, A, r)
points(sims, pch = 19, col = "green")






# 3D ####
A <- matrix(c(5,1,1, 1,3,1, 1,1,1), ncol = 3L)
r <- 2
# draw the ellipsoid
library(misc3d)
x <- seq(-1, 1, length.out = 50)
y <- seq(-1.5, 1.5, length.out = 50)
z <- seq(-2.7, 2.7, length.out = 50)
g <- as.matrix(expand.grid(x = x, y = y, z = z))
voxel <- 
  array(apply(g, 1L, function(v) t(v) %*% A %*% v), dim = c(50, 50, 50))
isosurface <- computeContour3d(voxel, max(voxel), r^2, x = x, y = y, z = z)
drawScene.rgl(makeTriangles(isosurface, alpha = 0.3))
# simulate and plot points on ellipsoid
library(rgl)
sims <- runif_on_ellipsoid(300, A, r)
points3d(sims)