gg <- function(x,r=1,p=2,usedeldir=TRUE)
{
	dx <- pdist(x,p=p)
   n <- nrow(dx)
   A <- matrix(0,nrow=n,ncol=n)
   if(is.vector(x)) x <- matrix(x,ncol=1)
   if(usedeldir && p==2 && ncol(x)==2 && 
	  (length(.find.package("deldir",quiet=TRUE))>0)){
	  require(deldir)
	  del <- deldir(x[,1],x[,2])
	  for(edge in 1:nrow(del$delsgs)){
	    i <- del$delsgs[edge,5]
	    j <- del$delsgs[edge,6]
		 d1 <- r*dx[i,j]/2
		 d <- pdist((x[i,]+x[j,])/2,x,d=ncol(x),p=p)
		 d[i] <- Inf
		 d[j] <- Inf
		 if(!any(d<d1)){
		    A[i,j] <- 1
		    A[j,i] <- 1
		 }
	  }
   }
   else{
	  for(i in 1:n){
		 for(j in setdiff(1:n,i)){
			d1 <- r*dx[i,j]/2
			d1 <- round(d1,10)
			d <- pdist((x[i,]+x[j,])/2,x,d=ncol(x),p=p)
			d <- round(d,10)
			d[i] <- Inf
			d[j] <- Inf
			if(!any(d<d1)){
				A[i,j] <- 1
				A[j,i] <- 1
			}
		 }
	  }
   }
	diag(A) <- 0
	g <- graph.adjacency(A,mode="undirected")
	g$r <- r
	g$layout <- x
	g$p <- p
   g
}

