\name{cccd}
\alias{cccd}
\alias{cccd.rw}
\alias{cccd.classify}
\alias{cccd.classifier}
\alias{cccd.classifier.rw}
\alias{cccd.multiclass.classifier}
\alias{cccd.multiclass.classify}
\alias{plot.cccd}
\alias{plot.cccdClassifier}
\title{
Class Cover Catch Digraph
}
\description{
Constructs a class cover catch digraph from points or interpoint distance matrices.
}
\usage{
cccd(x = NULL, y = NULL, dxx = NULL, dyx = NULL, method = NULL, 
     k = NA, algorithm = 'cover_tree')
cccd.rw(x=NULL,y=NULL,dxx=NULL,dyx=NULL,method=NULL,m=1,d=2)
cccd.classifier(x,y,dom.method='greedy',proportion=1,...)
cccd.classify(data, C,method=NULL)
cccd.classifier.rw(x,y,m=1,d=2)
cccd.multiclass.classifier(data, classes, dom.method='greedy',proportion=1,...)
cccd.multiclass.classify(data,C,method=NULL)
\method{plot}{cccd}(x, ..., plot.circles = FALSE, dominate.only = FALSE, 
          D = NULL, vertex.size = 2, vertex.label = NA, 
			 vertex.color = "SkyBlue2", dom.color = "Blue", 
			 ypch = 20, ycex = 1.5, ycol = 2, 
			 use.circle.radii = FALSE, balls = FALSE, 
			 ball.color = gray(0.8), square = FALSE, xlim, ylim)
\method{plot}{cccdClassifier}(x, ..., xcol=1,ycol=2,xpch=20,ypch=xpch,
                                balls=FALSE,add=FALSE)
}
\arguments{
    \item{x,y}{
	    the target class and non-target class points. Either x,y
		 or dxx,dyx must be provided. In the case of \code{plot}, x is
		 an object of class cccd.
    }
    \item{dxx,dyx}{
	    interpoint distances (x against x and y against x). If these
		 are not provided they are computed using x and y.
    }
  \item{method}{ the method used for the distance. 
     See \code{\link[proxy]{dist}}.}
  \item{dom.method,proportion}{ the method used for the domination set 
  computation, and the proportion of points required to dominate.  
  See \code{\link{dominate}}.}
  \item{k}{If given, \code{get.knn} is used from FNN to approximate
   the class cover catch graph. Each x covers no more than the \code{k} 
	nearest neighbors to it.
	This will be much faster and use less memory for large data sets, but 
	is an approximation unless \code{k} is sufficiently large.}
  \item{algorithm}{See \code{\link[FNN]{get.knn}}.}
	 \item{m}{slope of the null hypothesis curve}
	 \item{data}{data to be classified}
	 \item{d}{dimension of the data}
	 \item{classes}{class labels of the data}
	 \item{C}{cccd object}
  \item{plot.circles}{ logical. Plot the circles around the points if TRUE.}
  \item{dominate.only}{ logical. Only plot the digraph induced by the
      dominating set.}
  \item{D}{ a dominating set. Only used if dominate.only is TRUE. If
         dominate.only is TRUE and D is NULL, then \code{dominate}
			is called.}
  \item{vertex.size,vertex.color,vertex.label, dom.color}{parameters controling
		  the plotting of the vertices. \code{dom.color} is the color
		  of the vertices in the dominating set.}
  \item{balls, ball.color}{if \code{balls}=TRUE, the cover is plotted
        as filled balls, with
        \code{ball.color} controling their color. In the cass of
		  \code{cccdClassifier}, \code{balls} can be "x" or "y" indicating
		  that only one of the balls should be plotted.}
  \item{ypch,ycex,ycol}{ parameters for plotting the non-target points.}
  \item{xpch,xcol}{parameters for plotting the first class points.}
  \item{add}{logical. Should the classifier plot be added to an existing plot?}
  \item{use.circle.radii}{ logical. Ensure that the circles fit
      within the plot.}
  \item{square}{logical. Make the plot square.}
  \item{xlim,ylim}{if present, these control the plotting region.}
  \item{\dots}{ arguments passed to \code{cccd} or \code{plot}.}
}
\details{
   The class cover catch digraph is a graph with vertices defined by the
	points of \code{x} and edges defined according to the balls
	\eqn{B(x,d(x,Y))}. There is an edge between vertices
	\eqn{x_1,x_2} if \eqn{x_2\in B(x_1,d(x_1,Y))}. If \code{dyx} is not
	given and the method is 'euclidean', then \code{get.knnx} is used to
	find the nearest \code{y} to each x. If \code{k} is given, only
	the \code{k} nearest neighbors to each point are candidates for
	covering. Thus the cccd will be approximate, but the computation will
	(generally) be faster. Since \code{get.knn} uses Euclidean distance,
	these choices will only be valid for this distance metric.
	Since the graph will tend to be larger than
	otherwise, the dominating set computation will be slower, so one
	should trade-off speed of calculation, approximation, and the 
	\code{proportion} option to the dominating set (which can make that
	calculation faster at the cost of returning a subset of the dominating
	set).

}
\value{
   an object of class igraph. In addition, it contains the attributes:
	\item{R}{a vector of radii.}
	\item{Y}{the y vectors.}
	\item{layout}{the x vectors.}
	In the case of the classifier, the attributes are:
	\item{Rx,Ry}{vectors of radii.}
	\item{Cx,Cy}{the ball centers.}
 }
\note{
    The plotting assumes the cccd used Euclidean distance, and so the
	 balls/circles will be Euclidean balls/circles. If the 
	 method used in the distance was some other metric, you'll have
	 to plot the balls/circles yourself if you want them to be correct
	 on the plot.
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
\code{\link{ccd}}, \code{\link{rng}}, \code{\link{gg}}, \code{\link[proxy]{dist}},
\code{\link[FNN]{get.knn}}
\code{\link{dominate}}
}
\examples{
set.seed(456330)
z <- matrix(runif(1000),ncol=2)
ind <- which(z[,1]<.5 & z[,2]<.5)
x <- z[ind,]
y <- z[-ind,]
g <- cccd(x,y)
C <- cccd.classifier(x,y)
z2 <- matrix(runif(1000),ncol=2)
ind <- which(z2[,1]<.5 & z2[,2]<.5)
cls <- rep(0,nrow(z2))
cls[ind] <- 1
out <- cccd.classify(z2,C)
sum(out != cls)/nrow(z2)
\dontrun{
plot(g,plot.circles=TRUE,dominate.only=TRUE)
points(z2,col=2*(1-cls)+1,pch=20)
}
}
\keyword{graphs}
