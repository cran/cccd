\name{nng}
\alias{nng}
\title{ Nearest Neighbor Graphs}
\description{
  nearest neighbor, k-nearest neighbor, and mutual k-nearest neighbor 
  (di)graphs.
}
\usage{
nng(x = NULL, dx = NULL, k = 1, mutual = FALSE, method = NULL, 
    use.fnn = FALSE, algorithm = 'cover_tree')
}
\arguments{
  \item{x}{ a data matrix. Either x or dx is required}
  \item{dx}{ interpoint distance matrix}
  \item{k}{ number of neighbors}
  \item{mutual}{ logical. if true the neighbors must be mutual. See details.}
  \item{method}{ the method used for the distance. 
     See \code{\link[proxy]{dist}}}
  \item{use.fnn}{logical. If TRUE, \code{get.knn} from the FNN package
  is used to obtain the neighbors.}
  \item{algorithm}{see \code{\link[FNN]{get.knn}}.}
}
\details{
 a k-nearest neighbor graph is a digraph where each vertex is
 associated with an observation and there is a directed edge between
 the vertex and it's k nearest neighbors. A mutual k-nearest neighbor
 graph is a graph where there is an edge between x and y if x is
 one of the k nearest neighbors of y AND y is one of the k nearest
 neighbors of x.
}
\value{
  an object of class igraph with the extra attributes
  \item{layout}{ the x vectors.}
  \item{k,mutual,p}{arguments given to \code{nn}.}
}
\references{ 
D.J. Marchette, Random Graphs for Statistical Pattern Recognition,
John Wiley & Sons, 2004.
}
\author{ David J. Marchette david.marchette@navy.mil}
\seealso{
\code{\link[proxy]{dist}}
\code{\link[FNN]{get.knn}}
}

\examples{

x <- matrix(runif(100),ncol=2)

G1 <- nng(x,k=1)
\dontrun{
par(mfrow=c(2,2))
plot(G1)
}

G2 <- nng(x,k=2)
\dontrun{
plot(G2)
}

G5 <- nng(x,k=5)
\dontrun{
plot(G5)
}

G5m <- nng(x,k=5,mutual=TRUE)
\dontrun{
plot(G5m)
par(mfrow=c(1,1))
}

}
\keyword{ graphs }
