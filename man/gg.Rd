\name{gg}
\alias{gg}
\title{ Gabriel Graph}
\description{
  A Gabriel graph is one where the vertices are points and there
  is an edge between two points if the maximal ball between the points
  contains no other points.
}
\usage{
gg(x, r = 1, method = NULL, usedeldir = TRUE)
}
\arguments{
  \item{x}{ a matrix of observations.}
  \item{r}{ a multiplier on the ball radius.}
  \item{method}{ the method used for the distance. 
     See \code{\link[proxy]{dist}}}
  \item{usedeldir}{logical. Whether to use the deldir package or not.}
}
\details{
   places an edge between two points \eqn{i,j} if the ball centered between
	the points with radius \eqn{r d(i,j)/2} contains no other points.
}
\value{
  an object of class igraph. In addition it contains the attributes:
  \item{layout}{ the data.}
  \item{r,p}{arguments passed to \code{gg}}
}
\references{ 
K.R. Gabriel and R.R. Sokal,
A New Statistical Approach to Geographic Variation Analysis,
Systemic Zoology, 18, 259-278, 1969

D.J. Marchette, Random Graphs for Statistical Pattern Recognition,
John Wiley & Sons, 2004.

}
\author{ David J. Marchette}

\seealso{ \code{\link{rng}}, \code{\link[proxy]{dist}}}

\examples{
x <- matrix(runif(100),ncol=2)

g <- gg(x)
\dontrun{
plot(g)
}

}
\keyword{ graphs }
