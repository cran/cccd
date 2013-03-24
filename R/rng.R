rng <- function(x=NULL,dx=NULL,r=1, method=NULL,usedeldir=TRUE)
{
   if(is.null(dx)) {
	  if(is.null(x)) stop("One of x or dx must be given.")
	  dx <- as.matrix(dist(x,method=method))
   }
   n <- nrow(dx)
   A <- matrix(0,nrow=n,ncol=n)
   if(is.vector(x)) x <- matrix(x,ncol=1)
	if(!require(deldir,quietly=TRUE)) usedeldir <- FALSE
   if(usedeldir && ncol(x)==2){
	  del <- deldir(x[,1],x[,2])
	  for(edge in 1:nrow(del$delsgs)){
	     i <- del$delsgs[edge,5]
	     j <- del$delsgs[edge,6]
		 d <- min(apply(cbind(dx[i,-c(i,j)],dx[j,-c(i,j)]),1,max))
		 if(r*dx[i,j] < d){
		    A[i,j] <- 1
		    A[j,i] <- 1
		 }
	  }
   }
   else{
	  diag(dx) <- Inf
	  for(i in 1:n){
		 for(j in setdiff(1:n,i)){
			d <- min(apply(cbind(dx[i,-c(i,j)],dx[j,-c(i,j)]),1,max))
			if(r*dx[i,j] < d){
			   A[i,j] <- 1
			   A[j,i] <- 1
			}
		 }
	  }
   }
	diag(A) <- 0
	out <- graph.adjacency(A,mode="undirected")
	out$layout <- x
	out$r <- r
   out
}

