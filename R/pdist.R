
guess.dim <- function(x,y)
{
   if(is.vector(x)){
      if(is.null(y)){
	     d <- 1
	  }
	  else{
	     if(is.matrix(y)){
		    d <- ncol(y)
		 }
		 else{
		    if(length(x)==length(y)){
			   d <- length(x)
			}
			else{
			   d <- 1
			}
		 }
	  }
   }
   else{
      d <- ncol(x)
   }
   d
}

pdist <- function(x,y=NULL,d,p=2,w=NULL,arc=FALSE)
{
   if(missing(d)) d <- guess.dim(x,y)
   x <- matrix(x,ncol=d)
   nx <- nrow(x)
   if(!is.null(y)){
	  y <- matrix(y,ncol=d)
	  return(pdistxy(x,y,d,p,w))
   }
   out <- rep(0,nx^2)
   if(is.infinite(p)){
	  retval <- .C("pdistinf",
				   x = as.double(t(x)),
				   n = as.integer(nx),
				   d = as.integer(d),
				   w = as.double(w),
				   out=as.double(out),PACKAGE="cccd")
   }
   else if(!is.na(p)){
	  retval <- .C("pdist",
				   x = as.double(t(x)),
				   n = as.integer(nx),
				   d = as.integer(d),
				   p = as.double(p),
				   w = as.double(w),
				   out=as.double(out),PACKAGE="cccd")
   }
	else{
	   return(pdistcos(x,d=d,arc=arc))
	}
   return(matrix(retval$out,nrow=nrow(x)))
}

pdistcos <- function(x,d=ncol(x),arc=FALSE)
{
   x <- matrix(x,ncol=d)
   nx <- nrow(x)
   out <- rep(0,nx^2)
   retval <- .C("pdistcos",
				x = as.double(t(x)),
				n = as.integer(nx),
				d = as.integer(d),
				out=as.double(out),PACKAGE="cccd")
	if(arc) retval$out <- acos(1-retval$out)
   return(matrix(retval$out,nrow=nrow(x)))
}

pdistxy <- function(x,y,d,p=2,w=NULL)
{
   if(missing(d)){
	  d <- guess.dim(x,y)
   }
   x <- matrix(x,ncol=d)
   nx <- nrow(x)
   y <- matrix(y,ncol=d)
   ny <- nrow(y)

   out <- rep(0,nx*ny)
   if(is.infinite(p)){
	  retval <- .C("pdistxyinf",
				   x = as.double(t(x)),
				   y = as.double(t(y)),
				   nx = as.integer(nx),
				   ny = as.integer(ny),
				   d = as.integer(d),
				   w = as.double(w),
				   dist=as.double(out),PACKAGE="cccd")
   }
   else {
	  retval <- .C("pdistxy",
				   x = as.double(t(x)),
				   y = as.double(t(y)),
				   nx = as.integer(nx),
				   ny = as.integer(ny),
				   d = as.integer(d),
				   p = as.double(p),
				   w = as.double(w),
				   dist=as.double(out),PACKAGE="cccd")
   }
   matrix(retval$dist,byrow=TRUE,nrow=nx)
}

