cccd <- function(x=NULL,y=NULL,dxx=NULL,dyx=NULL,method=NULL)
{
   if(is.null(dxx) | is.null(dyx)){
		if(is.null(x) | is.null(y)) stop("either x,y or dxx,dyx must be given")
      dyx <- as.matrix(dist(y,x,method=method))
      dxx <- as.matrix(dist(x,method=method))
   }
   R <- apply(dyx,2,min)
	M <- matrix(as.integer(dxx<R),length(R))
	diag(M) <- 0
	g <- graph.adjacency(M,mode="directed")
	g$R <- R
	if(is.vector(x) || (dim(x)==1)) {
		x <- cbind(x,rep(0,length(x)))
		y <- cbind(y,rep(0,length(y)))
	}
	g$layout <- x
	g$Y <- y
	g$method <- method
	class(g) <- c("cccd",class(g))
	g
}

cccd.rw <- function(x=NULL,y=NULL,dxx=NULL,dyx=NULL,method=NULL,m=1,d=2)
{
   if(is.null(dxx) | is.null(dyx)){
		if(is.null(x) | is.null(y)) stop("either x,y or dxx,dyx must be given")
      dyx <- as.matrix(dist(y,x,method=method))
      dxx <- as.matrix(dist(x,method=method))
		d <- ncol(x)
	} 
   R <- rep(0,nrow(dxx))
	nx <- nrow(dxx)
	ny <- nrow(dyx)
	for(i in 1:nx){
		o <- order(c(dxx[i,],dyx[,i]))
	   rw <- cumsum(c(rep(1/nx,nx),rep(-1/ny,ny))[o])
		r <- sort(dxx[i,])
		R[i] <- r[which.max(rw[1:nx]-m*r^d)]
	}
	M <- matrix(as.integer(dxx<R),length(R))
	diag(M) <- 0
	g <- graph.adjacency(M,mode="directed")
	g$R <- R
	g$layout <- x
	g$Y <- y
	g$method <- method
	class(g) <- c("cccd",class(g))
	g
}

plot.cccd <- function(x,...,
                     plot.circles=FALSE,dominate.only=FALSE,D=NULL,
							vertex.size=2,vertex.label=NA,
							vertex.color="SkyBlue2",dom.color="Blue",
                     ypch=20,
							ycex=1.5,ycol=2,
							use.circle.radii=FALSE,
							balls=FALSE,
							ball.color=gray(.8),
							square=FALSE,
							xlim,ylim)

{
	g <- x
	class(g) <- "igraph"
	if(balls) plot.circles <- TRUE
	x <- g$layout
	n <- nrow(x)
	y <- g$Y
	if(is.null(y)){
		if(missing(xlim)){
			xlim <- range(x[,1])
		}
		if(missing(ylim)){
			ylim <- range(x[,2])
		}
	}
	else {
		if(missing(xlim)){
			xlim <- range(c(x[,1],y[,1]))
		}
		if(missing(ylim)){
			ylim <- range(c(x[,2],y[,2]))
		}
	}
	if(is.null(D)) D <- dominate(g)+1
	else D <- D+1
	col <- rep(vertex.color,n)
	col[D] <- dom.color
	vertex.color <- col
	r <- g$R
	if(use.circle.radii){
		 xlim <- range(c(xlim[2]+r,xlim[1]-r))
		 ylim <- range(c(ylim[2]+r,ylim[1]-r))
	}
	if(square){
	   xlim <- range(c(xlim,ylim))
	   ylim <- xlim
	}
	plot(g,xlim=xlim,ylim=ylim,vertex.size=vertex.size,rescale=FALSE,
	           vertex.label=vertex.label,
				  vertex.color=vertex.color,...)
	if(plot.circles){
		col <- rep(ifelse(balls,ball.color,vertex.color),n)
		if(dominate.only){
		   col[-D] <- NA
		}
		for(i in 1:nrow(x)){
			if(!is.na(col[i])){
				if(balls){
					draw.circle(x[i,1],x[i,2],r[i],border=vertex.color[i],col=col[i])
				} else {
					draw.circle(x[i,1],x[i,2],r[i],border=col[i])
				}
			}
		}
		if(balls){
			plot(g,xlim=xlim,ylim=ylim,vertex.size=vertex.size,rescale=FALSE,
						  vertex.label=vertex.label,
						  vertex.color=vertex.color,add=TRUE,...)
		}
	}
	if(!is.null(y)){
		points(y,pch=ypch,col=ycol,cex=ycex)
	}
}

cccd.classifier <- function(x,y,method=NULL)
{
	if(missing(y)){
	   if(is.list(x)){
		   y <- x$y
			x <- x$x
			if(is.null(x) | is.null(y))
			   stop("must provide either x and y or a list with attributes x and y")
		}
		else
			stop("must provide either x and y or a list with attributes x and y")
	}
   Gx <- cccd(x,y,method=method)
	Gy <- cccd(y,x,method=method)
	Dx <- dominate(Gx)
	Dy <- dominate(Gy)
	list(Rx=Gx$R[Dx],Ry=Gy$R[Dy],Cx=matrix(x[Dx,],ncol=ncol(x)),Cy=matrix(y[Dy,],ncol=ncol(y)))
}

cccd.classifier.rw <- function(x,y,m=1,d=2)
{
	if(missing(y)){
	   if(is.list(x)){
		   y <- x$y
			x <- x$x
			if(is.null(x) | is.null(y))
			   stop("must provide either x and y or a list with attributes x and y")
		}
		else
			stop("must provide either x and y or a list with attributes x and y")
	}
   Gx <- cccd.rw(x,y,m=m,d=d)
	Gy <- cccd.rw(y,x,m=m,d=d)
	Dx <- dominate(Gx)
	Dy <- dominate(Gy)
	list(Rx=Gx$R[Dx],Ry=Gy$R[Dy],Cx=matrix(x[Dx,],ncol=ncol(x)),Cy=matrix(y[Dy,],ncol=ncol(y)))
}

cccd.classify <- function(data,C,method=NULL)
{
	dx <- apply(t(t(as.matrix(dist(data,C$Cx,method=method)))/C$Rx),1,min)
	dy <- apply(t(t(as.matrix(dist(data,C$Cy,method=method)))/C$Ry),1,min)
	dx<dy
}

cccd.multiclass.classifier <- function(data,classes)
{
	cls <- unique(classes)
	nc <- length(cls)
	G <- list(0)
	D <- list(0)
	C <- list(0)
	R <- list(0)
	for(i in 1:nc){
	   z <- classes==cls[i]
		x <- data[z,]
		y <- data[!z,]
		G[[i]] <- cccd(x,y)
		D[[i]] <- dominate(G[[i]])
		C[[i]] <- matrix(x[D[[i]],],ncol=ncol(x))
		R[[i]] <- G[[i]]$R[D[[i]]]
	}
	list(G=G,D=D,C=C,R=R,classes=cls)
}

cccd.multiclass.classify <- function(data,C,method=NULL)
{
	nc <- length(C$R)
	if(is.vector(data)) data <- matrix(data,nrow=1)
	d <- matrix(0,nrow=nc,ncol=nrow(data))
	classes <- C$classes
	for(i in 1:nc){
		d[i,] <- apply(t(t(as.matrix(dist(data,C$C[[i]],method=method)))/C$R[[i]]),1,min)
	}
	z <- t(apply(d,2,function(x)x/sum(x)))
	list(probs=z,classes=classes[apply(z,1,which.min)])
}

