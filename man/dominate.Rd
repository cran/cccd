\name{dominate}
\alias{dominate}
\title{ Dominating Sets}
\description{
   find maximum dominating sets in (di)graphs.
}
\usage{
dominate(g, method = "greedy")
}
\arguments{
  \item{g}{ an adjacency matrix.}
  \item{method}{ one of "greedy","random","byRadius".}
}
\details{
   \code{dominate} is the main program which calls the others,
	as indicated by \code{method}. Greedy is the greedy dominating
	algorithm.  In the greedy method ties are broken by first index (a la 
	\code{which.max}).
	The byRadius method uses the radii to break ties while
	the random routine breaks ties randomly.
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
