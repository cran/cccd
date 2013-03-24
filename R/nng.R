
nng <- function(x=NULL,dx=NULL,k=1,mutual=FALSE,method=NULL)
{
   if(is.null(dx)) {
	  if(is.null(x)) stop("one of x or dx must be given")
	  dx <- as.matrix(dist(x,method=method))
   }
   n <- nrow(dx)
   A <- matrix(0,nrow=n,ncol=n)
   for(i in 1:n){
	  d <- sort(dx[i,])
      A[i,dx[i,]<=d[k+1]] <- 1
   }
   diag(A) <- 0
   if(mutual){
      for(i in 1:n){
	     A[i,] <- A[i,] & A[,i]
		 A[,i] <- A[i,]
	  }
   }
	if(mutual)
		out <- graph.adjacency(A,mode="undirected")
	else
		out <- graph.adjacency(A,mode="directed")
	out$k <- k
	out$mutual <- mutual
	out$layout <- x
   out
}

