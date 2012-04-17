\name{ccd}
\alias{ccd}
\alias{plot.ccd}
\title{ Cluster Catch Digraphs}
\description{
 construct the cluster catch digraph from a data matrix.
}
\usage{
ccd(data, m = 1, alpha = 0.05, sequential = TRUE, method = NULL)
\method{plot}{ccd}(x,...)
}
\arguments{
  \item{data}{ a matrix of observations.}
  \item{m}{ slope of the null hypothesis curve.}
  \item{alpha}{ alpha for the K-S test if \code{sequential=T}.}
  \item{sequential}{ use the sequential or non-sequential version.}
  \item{method}{ the method used for the distance. 
     See \code{\link[proxy]{dist}}.}
  \item{x}{an object of class ccd.}
  \item{\dots}{arguments passed to \code{plot.cccd}.}
}
\details{
  cluster cover digraph. \code{plot.ccd} is just a call to \code{plot.cccd}.
}
\value{
 an object of class igraph. In addition, this contains the attributes:
 \item{R}{the radii.}
 \item{stats}{ the K-S statistics.}
 \item{layout}{the data vectors.}
 \item{walks}{the y-values of the random walks.}
 \item{fs}{the null hypothesis curve.}
 \item{A}{ the adjacency matrix.}
 \item{m,alpha}{arguments passed to \code{ccd}.}
}
\references{ 
D.J. Marchette, Random Graphs for Statistical Pattern Recognition,
John Wiley & Sons, 2004.

}
\author{ David J. Marchette david.marchette@navy.mil}

\seealso{ \code{\link{cccd}}
}

\examples{

x <- matrix(rnorm(100),ncol=2)
G <- ccd(x)
\dontrun{
plot(G)
}
}
\keyword{ graphs }
