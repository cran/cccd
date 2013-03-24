\name{prune}
\alias{prune}
\title{ Prune Points}
\description{
   a nearest neighbor pruning using neighborhood graphs.
}
\usage{
prune(x, classes, prox = "Gabriel", ignore.ties = TRUE, ...)
}
\arguments{
  \item{x}{ a data matrix.}
  \item{classes}{ a vector of class labels.}
  \item{prox}{ type of proximity graph.}
  \item{ignore.ties}{ do not prune if there is a tie vote.}
  \item{\dots}{ arguments passed to the proximity graph.}
}
\details{
   First a proximity graph is computed on the data. Then points
	are marked if their neighbors have a different class than they
	do: if the most common class among the neighbors is different
	than the point. Then all marked points are removed.
}
\value{
  A list with attributes:
  \item{x}{the pruned data.}
  \item{v}{the indices of the retained data.}
  \item{g}{the proximity graph.}
}
\references{ \url{http://www.bic.mni.mcgill.ca/users/crisco/pgedit/}}
\author{ David J. Marchette, david.marchette@navy.mil}

\keyword{ graphs }
