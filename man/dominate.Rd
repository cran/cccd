\name{dominate}
\alias{dominate}
\alias{dominate.greedy}
\alias{dominate.byR}
\alias{dominate.systematic}
\alias{dominate.random.sample}
\title{ Dominating Sets}
\description{
   find maximum dominating sets in (di)graphs.
}
\usage{
dominate(g, method = "greedy", verbose = TRUE)
dominate.greedy(g, weight=NULL)
dominate.byR(g)
dominate.random.sample(g)
}
\arguments{
  \item{g}{ an adjacency matrix.}
  \item{method}{ one of "greedy","systematic","random","byR".}
  \item{weight}{ weight vector.}
  \item{verbose}{ logical. In the case of systematic, whether to print out
                  information as it goes.}
}
\details{
   \code{dominate} is the main program which calls the others,
	as indicated by \code{method}. Greedy is the greedy dominating
	algorithm. The weight vector is used to break ties: in the case
	of a tie the vertex with the largest weight is chosen. The
	"by radius" algorithm uses the vector \code{R} as the criterion
	in the greedy algorithm, rather than degree. 
	The random routine simply uses a random weighting instead of weighting
	by \code{R}. The systematic routine computes a true minimum dominating
	set, but may take a long time as it is not optimized and uses brute force.
	It first computes a greedy dominating set to reduce the search, but
	that's all.
}
\value{
  a vector of vertices corresponding to the dominating set.
  Note: just like the vertex labels in igraph, these are 0-based.
}
\references{ 
T.W. Haynes, S.T. Hedetniemi and P.J. Slater,
Fundamentals of Domination in Graphs,
Marcel Dekker,
1998,
}
\author{ David J. Marchette david.marchette@navy.mil}

\examples{

x <- matrix(runif(100),ncol=2)
y <- matrix(runif(100,-2,2),ncol=2)
G <- cccd(x,y)
D <- dominate(G)
\dontrun{
plot(G,balls=TRUE,D=D)
}

}
\keyword{ graphs }
