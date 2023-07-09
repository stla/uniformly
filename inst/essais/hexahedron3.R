makeHexahedron <- function(p0, p1, p2, p3, p4, p5, p6, p7) {
  cbind(p0, p1, p2, p3, p4, p5, p6, p7)
}

hexahedron <- makeHexahedron(
  p0 = c(1.7, 1.7, 0),
  p1 = c(2, 0, 0),
  p2 = c(2, 2, 0),
  p3 = c(0, 2, 0),
  p4 = c(0, 2, 2),
  p5 = c(0, 0, 2),
  p6 = c(2, 0, 2),
  p7 = c(2, 2, 2)
)

library(rgl)
plotHexahedron <- function(hexahedron, alpha = 1) {
  faces <- cbind(
    c(4L, 3L, 2L, 1L),
    c(5L, 6L, 7L, 8L),
    c(1L, 2L, 7L, 6L),
    c(6L, 5L, 4L, 1L),
    c(2L, 3L, 8L, 7L),
    c(3L, 4L, 5L, 8L)
  )
  mesh <- qmesh3d(
    vertices = hexahedron,
    indices  = faces
  )
  shade3d(mesh, color = "yellow", alpha = alpha)
  wire3d(mesh, color = "black")
}

library(abind)

tetrahedra_hexahedron <- function(hexahedron) {
  vs1 <- hexahedron[, c(1L, 5L, 3L, 7L)]
  vs2 <- hexahedron[, c(7L, 1L, 2L, 3L)]
  vs3 <- hexahedron[, c(6L, 1L, 7L, 5L)]
  vs4 <- hexahedron[, c(8L, 7L, 3L, 5L)]
  vs5 <- hexahedron[, c(4L, 1L, 3L, 5L)]
  abind(vs1, vs2, vs3, vs4, vs5, along = 3L)  
}

volume_hexahedron <- function(hexahedron) {
  vs <- tetrahedra_hexahedron(hexahedron)
  vs1 <- vs[, , 1L]
  vs2 <- vs[, , 2L]
  vs3 <- vs[, , 3L]
  vs4 <- vs[, , 4L]
  vs5 <- vs[, , 5L]
  vol1 <- volume_tetrahedron(vs1[, 1L], vs1[, 2L], vs1[, 3L], vs1[, 4L])
  vol2 <- volume_tetrahedron(vs2[, 1L], vs2[, 2L], vs2[, 3L], vs2[, 4L])
  vol3 <- volume_tetrahedron(vs3[, 1L], vs3[, 2L], vs3[, 3L], vs3[, 4L])
  vol4 <- volume_tetrahedron(vs4[, 1L], vs4[, 2L], vs4[, 3L], vs4[, 4L])
  vol5 <- volume_tetrahedron(vs5[, 1L], vs5[, 2L], vs5[, 3L], vs5[, 4L])
  vol1 + vol2 + vol3 + vol4 + vol5  
}

runif_in_hexahedron <- function(n, hexahedron) {
  # tetrahedra ####
  vs <- tetrahedra_hexahedron(hexahedron)
  vs1 <- vs[, , 1L]
  vs2 <- vs[, , 2L]
  vs3 <- vs[, , 3L]
  vs4 <- vs[, , 4L]
  vs5 <- vs[, , 5L]
  vol1 <- volume_tetrahedron(vs1[, 1L], vs1[, 2L], vs1[, 3L], vs1[, 4L])
  vol2 <- volume_tetrahedron(vs2[, 1L], vs2[, 2L], vs2[, 3L], vs2[, 4L])
  vol3 <- volume_tetrahedron(vs3[, 1L], vs3[, 2L], vs3[, 3L], vs3[, 4L])
  vol4 <- volume_tetrahedron(vs4[, 1L], vs4[, 2L], vs4[, 3L], vs4[, 4L])
  vol5 <- volume_tetrahedron(vs5[, 1L], vs5[, 2L], vs5[, 3L], vs5[, 4L])
  volumes <- c(vol1, vol2, vol3, vol4, vol5)
  probs <- volumes / sum(volumes)
  # sampling
  indices <- sample.int(5L, size = n, replace = TRUE, prob = probs)
  sims <- matrix(NA_real_, nrow = n, ncol = 3L)
  for(i in 1L:n) {
    th <- vs[, , indices[i]]
    sims[i, ] <- 
      runif_in_tetrahedron(1L, th[, 1L], th[, 2L], th[, 3L], th[, 4L])
  }
  sims
}

sims <- runif_in_hexahedron(5000, hexahedron)
plotHexahedron(hexahedron, alpha = 0.2)
points3d(sims)

in_tetrahedron <- function(M, tetrahedron) {
  A <- rbind(tetrahedron, 1)
  u <- solve(A, c(M, 1))
  all(u >= 0)
}
in_hexahedron <- function(M, hexahedron) {
  vs <- tetrahedra_hexahedron(hexahedron)
  th1 <- vs[, , 1L]
  th2 <- vs[, , 2L]
  th3 <- vs[, , 3L]
  th4 <- vs[, , 4L]
  th5 <- vs[, , 5L]
  in_tetrahedron(M, th1) || in_tetrahedron(M, th2) || in_tetrahedron(M, th3) || 
    in_tetrahedron(M, th4) || in_tetrahedron(M, th5)
}