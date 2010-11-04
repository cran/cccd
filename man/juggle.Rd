\name{juggling}
\alias{juggle}
\alias{juggle.classify}
\title{
Juggling
}
\description{
a resampled version of the CCCD classifier.
}
\usage{
juggle(data, classes, sampled = TRUE, sample.dim = FALSE, 
       num = 100, sample.proportion = 0.1, k = 2)
juggle.classify(data,J,tdata,indices)
}
\arguments{
  \item{data,tdata}{
training data from which to build the classifier. In the case
of \code{juggle.classify}, \code{tdata} is the training data and
\code{data} is the test data.
}
  \item{classes}{
class labels.
}
  \item{sampled}{
whether the data are subsampled.
}
  \item{sample.dim}{
if TRUE, the dimensions (variates) are also sampled.
}
  \item{num}{
number of juggles (resamples).
}
  \item{sample.proportion}{
proportion of the data to sample. If 1 or greater, the data are sampled
with replacement.
}
  \item{k}{
number of variates to sample when \code{sample.dim} is TRUE.
}
  \item{J}{
the juggled classifier.
}
  \item{indices}{
the indices of the juggles to use.
}
}
\details{
The idea of juggling is to sample the data, compute a CCCD classifier,
then repeat. The resampling is controled by the two sampling variables,
which basically determine whether the data are sampled with replacement,
or whether a subsample is used. If \code{sample.dim} is TRUE, the 
variates are also sampled, with \code{k} indicating how many are sampled.
}
\value{
\code{juggle.classify} returns a matrix holding the classification 
probabilities for each observation in \code{data}.
a list consisting of:
\item{S }{the dominating sets.}
\item{R }{the radii.}
\item{dimension }{the dimension of the data.}
\item{vars}{in the case of \code{sample.dim}=TRUE, the variables sampled
            each time.}
Only the indicies into the training data are stored in \code{J}, which
is why the classifier requires the original training data in \code{tdata}.
}
\author{
David J. Marchette, david.marchette@navy.mil
}

\seealso{
\code{\link{cccd}}
}
\keyword{ math }
\keyword{ multivariate }
