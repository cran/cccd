\name{cccd-internal}
\alias{cccd-internal}
\alias{ccd.nonsequential}
\alias{ccd.sequential}
\alias{guess.dim}
\alias{pdistxy}
\alias{vote}
\alias{juggle1}
\alias{juggle2}
\alias{juggle1.classify}
\alias{juggle2.classify}
\title{
Internal functions for graph package
}
\description{
Internal functions that should not be called by the user.
}
\usage{
ccd.sequential(data,m=1,alpha=0.05)
ccd.nonsequential(data,m=1)
guess.dim(x,y)
vote(x)
pdistxy(x,y,d,p=2,w=NULL)
juggle1(data,classes,sampled=TRUE,
                   num=100,sample.proportion=0.1)
juggle2(data,classes,sampled=TRUE,
                   num=100,sample.proportion=0.1,
						 k=2)
juggle1.classify(data,J,tdata,indices)
juggle2.classify(data,J,tdata,indices)
}
\arguments{
    \item{data}{
	   a matrix of observations.
    }
    \item{m,alpha}{
	    arguments to \code{ccd}.
    }
	\item{x,y}{
	   vectors or matrices.
	}
	\item{d,p,w}{arguments to \code{pdist}}
	\item{sampled,classes,num,sample.proportion,k,J,tdata,indices}{
	      see \code{juggle}
	}
}
\details{
\itemize{
\item{ccd.sequential,ccd.nonsequential}{ sequential and nonsequential
versions of \code{ccd}. Called by \code{ccd}.}
\item{guess.dim}{ tries to guess the dimension for \code{pdist}.}
\item{vote}{called by \code{prune}. \code{pdistxy} is called by
\code{pdist} when computing the distances between two data sets.}
}
}
\author{
David J. Marchette, david.marchette@navy.mil
}
\seealso{
\code{\link{ccd}},\code{\link{pdist}}, 
      \code{\link{juggle}}, \code{\link{prune}},
}
\keyword{math}
