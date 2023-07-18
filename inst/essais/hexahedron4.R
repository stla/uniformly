library(uniformly)

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

sims <- runif_in_hexahedron(6000, hexahedron)

in_tetrahedron <- function(M, tetrahedron) {
  A <- rbind(tetrahedron, 1)
  u <- solve(A, c(M, 1))
  all(u >= 0)
}
in_hexahedron <- function(M, hexahedron) {
  vs <- uniformly:::tetrahedra_hexahedron(hexahedron)
  th1 <- vs[, , 1L]
  th2 <- vs[, , 2L]
  th3 <- vs[, , 3L]
  th4 <- vs[, , 4L]
  th5 <- vs[, , 5L]
  in_tetrahedron(M, th1) || in_tetrahedron(M, th2) || in_tetrahedron(M, th3) || 
    in_tetrahedron(M, th4) || in_tetrahedron(M, th5)
}

sims_cube <- runif_in_cube(1000, d = 3, O = c(1, 1, 1), r = 1)
inhx <- apply(sims_cube, 1L, function(M) in_hexahedron(M, hexahedron))
8 * mean(inhx)
volume_hexahedron(hexahedron)
