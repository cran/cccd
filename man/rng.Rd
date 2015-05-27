\name{rng}
\alias{rng}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Relative Neighborhood Graph.}
\description{
   the relative neighborhood graph defined by a set of points.
}
\usage{
rng(x=NULL, dx=NULL, r = 1, method = NULL, usedeldir = TRUE, open = TRUE, k = NA,
    algorithm = 'cover_tree')
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ a data matrix. Either \code{x} or \code{dx} must be provided.}
  \item{dx}{ an interpoint distance matrix.}
  \item{r}{ a multiplier to grow the balls.}
  \item{method}{ the method used for the distance. 
     See \code{\link[proxy]{dist}}}
  \item{usedeldir}{a logical. If true and the data are two dimensional
        and the deldir package is installed, the Delaunay triangularization
		  is first computed, and this is used to compute the relative
		  neighborhood graph.}
  \item{open}{logical. If TRUE, open balls are used in the definition.}
  \item{k}{If given, \code{get.knn} is used from FNN to approximate
   the relative neighborhood graph. Only the \code{k} nearest neighbors to
	the points are used to determine whether an edge should be made or not.
	This will be much faster and use less memory for large data sets, but 
	is an approximation unless \code{k} is sufficiently large.}
  \item{algorithm}{See \code{\link[FNN]{get.knn}}.}
}
\details{
   the relative neighborhood graph is defined in terms of balls
	centered at observations. For two observations, the balls are
	set to have radius equal to the distance between the observations
	(or \code{r} times this distance if \code{r} is not 1). There is
	an edge between the vertices associated with the observations if 
	and only if there are no vertices in the lune defined by the
	intersection of the balls.

	The flag \code{open} should make no difference for most applications,
	but there are very specific cases (see the example section below)
	where setting it to be TRUE will give the wrong answer (thanks to
	Luke Mathieson for pointing this out to me).
}
\value{
	an object of class igraph, with the additional attributes
   \item{layout}{the x matrix.}
	\item{r,p}{arguments given to \code{rng}.}
}
\references{ 

J.W. Jaromczyk and G.T. Toussaint,
"Relative neighborhood graphs and their relatives",
Proceedings of the IEEE, 
80, 1502-1517, 1992.

G.T. Toussaint,
"A Graph-Theoretic Primal Sketch",
Computational Morphology, 229-260, 1988.

D.J. Marchette, Random Graphs for Statistical Pattern Recognition,
John Wiley & Sons, 2004.

}
\author{ David J. Marchette david.marchette@navy.mil}

\seealso{\code{\link{gg}},\code{\link{cccd}},\code{\link{ccd}},
\code{\link[proxy]{dist}}
\code{\link[FNN]{get.knn}}
}

\examples{
x <- matrix(runif(100),ncol=2)

g <- rng(x)
\dontrun{
plot(g)
}

## Example using 'open':
g <- graph.full(5,directed=FALSE)

g1 <- rng(x=get.adjacency(g,sparse=FALSE),open=TRUE)
ecount(g1)
g2 <- rng(x=get.adjacency(g,sparse=FALSE),open=FALSE)
graph.isomorphic(g2,g)


}
\keyword{ graphs }
