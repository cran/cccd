\name{plotCCCD}
\alias{plotCCCD}
\alias{plotCCD}
\title{ Plot a class cover catch digraph.}
\description{
   plot routine for class cover catch digraphs.
}
\usage{
plotCCCD(g, plot.circles = FALSE, dominate.only = FALSE, 
          D = NULL, vertex.size = 2, vertex.label = NA, 
			 vertex.color = "SkyBlue2", dom.color = "Blue", 
			 ypch = 20, ycex = 1.5, ycol = 2, 
			 use.circle.radii = FALSE, balls = FALSE, 
			 ball.color = gray(0.8), square = FALSE, xlim, ylim, ...)
plotCCD(g, ...)
}
\arguments{
  \item{g}{ a CCCD graph or CCD clustering}
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
        \code{ball.color} controling their color}.
  \item{ypch,ycex,ycol}{ parameters for plotting the non-target points.}
  \item{use.circle.radii}{ logical. Ensure that the circles fit
      within the plot.}
  \item{square}{logical. Make the plot square.}
  \item{xlim,ylim}{if present, these control the plotting region.}
  \item{\dots}{ arguments passed to \code{plot} or \code{plotCCCD}.}
}
\details{
   this plots the class cover catch digraph, with the target and
	non-target points, or a cluster catch digraph. \code{plotCCD}
	is just a call to \code{plotCCCD}.
}
\author{ David J. Marchette david.marchette@navy.mil}

\keyword{ multivariate }
