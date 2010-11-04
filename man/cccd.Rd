\name{cccd}
\alias{cccd}
\alias{cccd.rw}
\alias{cccd.classify}
\alias{cccd.classifier}
\alias{cccd.classifier.rw}
\alias{cccd.multiclass.classifier}
\alias{cccd.multiclass.classify}
\title{
Class Cover Catch Digraph
}
\description{
Constructs a class cover catch digraph from points or interpoint distance matrices.
}
\usage{
cccd(x = NULL, y = NULL, dxx = NULL, dyx = NULL, p = 2, d = 2)
cccd.rw(x=NULL,y=NULL,dxx=NULL,dyx=NULL,p=2,d=2,m=1)
cccd.classifier(x,y)
cccd.classify(data, C)
cccd.classifier.rw(x,y,m=1,d=2)
cccd.multiclass.classifier(data, classes)
cccd.multiclass.classify(data,C)
}
\arguments{
    \item{x,y}{
	    the target class and non-target class points. Either x,y
		 or dxx,dyx must be provided.
    }
    \item{dxx,dyx}{
	    interpoint distances (x against x and y against x). If these
		 are not provided they are computed using x and y.
    }
	 \item{p}{
	    p in the \eqn{L_p} distance, if \code{dxx} or \code{dyx} are missing.}
	 \item{d}{dimension of the data.}
	 \item{m}{slope of the null hypothesis curve}
	 \item{data}{data to be classified}
	 \item{classes}{class labels of the data}
	 \item{C}{cccd object}
}
\details{
   The class cover catch digraph is a graph with vertices defined by the
	points of \code{x} and edges defined according to the balls
	\eqn{B(x,d(x,Y))}. There is an edge between vertices
	\eqn{x_1,x_2} if \eqn{x_2\in B(x_1,d(x_1,Y))}.
}
\value{
   an object of class igraph. In addition, it contains the attributes:
	\item{R}{a vector of radii.}
	\item{Y}{the y vectors.}
	\item{layout}{the x vectors.}
 }
\references{
D.J. Marchette,
"Class Cover Catch Digraphs",
Wiley Interdisciplinary Reviews: Computational Statistics, 
2, 171-177, 2010.

D.J. Marchette, Random Graphs for Statistical Pattern Recognition,
John Wiley & Sons, 2004.

C.E. Priebe, D.J. Marchette, J. DeVinney and D. Socolinsky, 
"Classification Using Class Cover Catch Digraphs",
Journal of Classification,
20, 3-23, 2003.
}
\author{
David J. Marchette, david.marchette@navy.mil
}
\seealso{
\code{\link{ccd}}, \code{\link{rng}}, \code{\link{gg}}
}
\keyword{math}
\keyword{multivariate}
