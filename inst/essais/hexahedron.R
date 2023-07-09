p0 <- c(0, 0, 0)
p1 <- c(2, 0, 0)
p2 <- c(2, 2, 0)
p3 <- c(0, 2, 0)
p4 <- c(0, 2, 2)
p5 <- c(0, 0, 2)
p6 <- c(2, 0, 2)
p7 <- c(2, 2, 2)

vs <- cbind(p0, p1, p2, p3, p4, p5, p6, p7)
faces <- cbind(
  c(4L, 3L, 2L, 1L),
  c(5L, 6L, 7L, 8L),
  c(1L, 2L, 7L, 6L),
  c(6L, 5L, 4L, 1L),
  c(2L, 3L, 8L, 7L),
  c(3L, 4L, 5L, 8L)
)

library(rgl)
mesh <- qmesh3d(
  vertices = vs,
  indices  = faces
)

shade3d(mesh, color = "yellow")
wire3d(mesh, color = "black")


# tetrahedra ####
faces <- cbind(
  c(1L, 2L, 3L),
  c(3L, 2L, 4L),
  c(4L, 2L, 1L),
  c(1L, 3L, 4L)
)

pts <- vs
vs1 <- pts[, c(1L, 5L, 3L, 7L)]
vs2 <- pts[, c(7L, 1L, 2L, 3L)]
vs3 <- pts[, c(6L, 1L, 7L, 5L)]
vs4 <- pts[, c(8L, 7L, 3L, 5L)]
vs5 <- pts[, c(4L, 1L, 3L, 5L)]

tth1 <- tmesh3d(
  vertices = vs1,
  indices  = faces
)
tth2 <- tmesh3d(
  vertices = vs2,
  indices  = faces
)
tth3 <- tmesh3d(
  vertices = vs3,
  indices  = faces
)
tth4 <- tmesh3d(
  vertices = vs4,
  indices  = faces
)
tth5 <- tmesh3d(
  vertices = vs5,
  indices  = faces
)

open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
shade3d(tth1, color = "red")
shade3d(tth2, color = "green")
shade3d(tth3, color = "blue")
shade3d(tth4, color = "yellow")
shade3d(tth5, color = "gray")
wire3d(mesh, color = "black")

# sampling ####
library(uniformly)
vol1 <- volume_tetrahedron(vs1[, 1L], vs1[, 2L], vs1[, 3L], vs1[, 4L])
vol2 <- volume_tetrahedron(vs2[, 1L], vs2[, 2L], vs2[, 3L], vs2[, 4L])
vol3 <- volume_tetrahedron(vs3[, 1L], vs3[, 2L], vs3[, 3L], vs3[, 4L])
vol4 <- volume_tetrahedron(vs4[, 1L], vs4[, 2L], vs4[, 3L], vs4[, 4L])
vol5 <- volume_tetrahedron(vs5[, 1L], vs5[, 2L], vs5[, 3L], vs5[, 4L])
volumes <- c(vol1, vol2, vol3, vol4, vol5)
volTotal <- sum(volumes)
probs <- volumes / volTotal

library(abind)
tetrahedra <- abind(vs1, vs2, vs3, vs4, vs5, along = 3L)

nsims <- 1000L
indices <- sample.int(5L, size = nsims, replace = TRUE, prob = probs)
sims <- matrix(NA_real_, nrow = nsims, ncol = 3L)
for(i in 1L:nsims) {
  th <- tetrahedra[, , indices[i]]
  sims[i, ] <- runif_in_tetrahedron(1L, th[, 1L], th[, 2L], th[, 3L], th[, 4L])
}

points3d(sims)
