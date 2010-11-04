\name{pdist}
\alias{pdist}
\alias{pdistcos}
\title{
Lp distances between vectors
}
\description{
computes the Lp distance between vectors. When called with a single
argument this is similar to \code{\link{dist}}. When called with
two arguments, it returns a matrix containing the interpoint distances
between the rows of the matrices.
}
\usage{
pdist(x,y=NULL,d,p=2,w=NULL,arc=FALSE)
pdistcos(x,d=ncol(x),arc=FALSE)
}
\arguments{
    \item{x}{
	   a vector or matrix.
    }
    \item{y}{
	   a vector or matrix. If y is NULL, pdist returns dist(x).
	   Otherwise it returns the distances between x and y.
    }
    \item{d}{
	   the dimension of the data (if missing, it guesses it from x and y).
    }
    \item{p}{
	   the p in Lp. p=Inf returns \eqn{L_\infty}. p=NA will call \code{pdistcos}.
    }
    \item{w}{
	   optional weight vector for weighted distances.
    }
    \item{arc}{logical. Should the arc-length be returned instead of the
	 cosine distance?
	 }
}
\keyword{utilities}
\details{
   This computes the \eqn{L_p} distance between
   the points in \code{x}, or, if \code{y} is not NULL, between pairs
   of \code{x} and \code{y}. If \code{p} is \code{Inf}, the \eqn{L_{\infty}}
	distance is computed. If \code{w} is given, it must be a vector of length
	equal to \code{d}, and is used to weight the contribution of the variates
	to the distance calculation. \code{pdistcos} returns the cosine
	dissimilarity, which is the normalized dot product of the vectors. If
	\code{arc}=TRUE, the arccos is returned, which is the arc-length distance.
	If \code{p}=NA, \code{pdist} will return \code{pdistcos} on the data.
}
\value{
   a matrix of distances. 
}
\seealso{
\code{\link{dist}}
}
\author{
David J. Marchette, david.marchette@navy.mil
}
\examples{
x <- runif(10)

a <- pdist(x)

# a is a 10x10 matrix of distances

y <- runif(20)
a <- pdist(x,y)

# a is now a 10x20 matrix of distances (it infered that d=1)

z <- runif(10)
b <- pdist(x,z)

# b is a 1x1 matrix (it guessed that d=10)

b <- pdist(x,z,d=1)

# now b is a 10x10 matrix 

z <- matrix(runif(20),ncol=10)

b <- pdist(x,z)

# now b is 1x2
}
