library(PlaneGeometry)
library(uniformly)

a <- 5; b <- 2

ell <- Ellipse$new(c(0, 0), a, b, 0)
perimeter <- ell$perimeter()


nsims0 <- 1e7
ssims <- runif_on_sphere(nsims0, 2)
allsims <- cbind(a*ssims[,1], b*ssims[,2])
p <- sqrt((a*ssims[,2])^2 + (b*ssims[,1])^2) / a
e <- runif(nsims0) < p 
sims <- allsims[e, ]
summary(sims)


alpha1 <- 10 * pi/180
alpha2 <- 60 * pi/180

mean(atan2(sims[,2], sims[,1]) > alpha1 & atan2(sims[,2], sims[,1]) < alpha2) * perimeter

ellArc <- EllipticalArc$new(ell, alpha1=alpha1, alpha2=alpha2, degrees = FALSE)
arc <- ellArc$length()
arc


sims2 <- runif_on_ellipsoid(100000, diag(c(1/a^2,1/b^2)), r = 1)
summary(sims2)
mean(atan2(sims2[,2], sims2[,1]) > alpha1 & atan2(sims2[,2], sims2[,1]) < alpha2) * perimeter

plot(ecdf(sims[, 1]))
lines(ecdf(sims2[, 1]), col = "red")

plot(ecdf(sims[, 2]))
lines(ecdf(sims2[, 2]), col = "blue")

cor(sims)
cor(sims2)
