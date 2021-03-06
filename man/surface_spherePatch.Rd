% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surfaces.R
\name{surface_spherePatch}
\alias{surface_spherePatch}
\title{Sphere patch surface}
\usage{
surface_spherePatch(r, phi1, phi2, theta1, theta2)
}
\arguments{
\item{r}{radius}

\item{phi1, phi2}{numbers defining the latitudinal angle range}

\item{theta1, theta2}{numbers defining the longitudinal angle range}
}
\value{
The surface of the sphere patch.
}
\description{
Surface of a sphere patch.
}
\details{
A sphere patch is the part of the sphere whose polar angles
\code{theta} and \code{phi} satisfy
\code{0 <= theta1 <= theta <= theta2 <= 2*pi} and
\code{0 <= phi1 <= phi <= phi2 <= pi}.
}
\examples{
# surface of the first orthant:
surface_spherePatch(r=1, phi1=0, phi2=pi/2, theta1=0, theta2=pi/2)
surface_stri(r=1, c(1,0,0), c(0,1,0), c(0,0,1))
}
\seealso{
\code{\link{surface_stri}} for the surface of a spherical triangle.
}
