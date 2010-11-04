
# This software developed by David Marchette
# Naval Surface Warfare Center, Dahlgren
# Division, Code B10.  It may be freely used and 
# modified. All warranties, express or implied, are
# disclaimed. 

cccd <- function(dxx,dyx,beta=0,verbose=F,shrink=T)
{

   require(utilities)

   n1 <- nrow(dxx)
   n2 <- nrow(dyx)
   if(beta>n2) beta <- n2

   retval <- .C("cccd",
		dx = as.double(dxx),
		dy = as.double(dyx),
		nx = as.integer(n1),
		ny = as.integer(n2),
		R = double(n1),
		C = integer(n1*n1),
		BETA = as.integer(beta),
		SHRINK = as.integer(shrink),
		VERBOSE = as.integer(verbose))

	  R <- retval$R

	  C <- matrix(retval$C,nrow=n1,byrow=T)

   out <- list(Indices=1:n1,Radii=R,AdjacencyMatrix=C,NumCovered=apply(C,1,sum),Covered=C)
   class(out) <- "cccd"
   out
}

dominate.cccd <- function(G,alpha=0,verbose=F)
{
   C <- G$AdjacencyMatrix
   N <- G$Num
   A <- NULL
   NC <- 0
   nx <- nrow(C)
   if(alpha>nx) alpha<-nx
   num <- nx-alpha

   covered <- rep(0,nx)

   while(NC<num){
      i <- biggestplace(N)
	  A <- c(A,i)
	  covered <- covered | C[i,]
	  NC <- sum(covered)
	  Q <- C[i,]
	  C <- t(apply(C,1,function(x,y){x[y>0] <- 0; x},Q))
	  N <- apply(C,1,sum)
	  i <- i+1
   }
   nw <- length(A)
   AM <- G$AdjacencyMatrix[A,A]
   if(is.matrix(AM)) NC <- apply(G$AdjacencyMatrix[A,],1,sum)
   else NC <- sum(G$AdjacencyMatrix)
   out <- list(Indices=A,Radii=G$Radii[A],AdjacencyMatrix=AM,
               NumCovered=NC,
               Covered=G$Covered[A,])
   class(out) <- "cccd"
   out
}

dominate.depth.cccd <- function(G,dxx,alpha=0,verbose=F)
{
   C <- G$AdjacencyMatrix
   N <- G$Num
   A <- NULL
   NC <- 0
   nx <- nrow(C)
   if(alpha>nx) alpha<-nx
   num <- nx-alpha

   dep <- 1/(1+apply(dxx,1,sum))

   covered <- rep(0,nx)

   while(NC<num){
      i <- biggestplace(N)
	  w <- N==N[i]
	  # break ties by largest depth
	  if(sum(w)>1){
		 qdep <- dep
		 qdep[!w] <- -1
		 i <- biggestplace(qdep)
	  }
	  A <- c(A,i)
	  covered <- covered | C[i,]
	  NC <- sum(covered)
	  Q <- C[i,]
	  C <- t(apply(C,1,function(x,y){x[y>0] <- 0; x},Q))
	  N <- apply(C,1,sum)
	  i <- i+1
   }
   nw <- length(A)
   AM <- G$AdjacencyMatrix[A,A]
   if(is.matrix(AM)) NC <- apply(G$AdjacencyMatrix[A,],1,sum)
   else NC <- sum(G$AdjacencyMatrix)
   out <- list(Indices=A,Radii=G$Radii[A],AdjacencyMatrix=AM,
               NumCovered=NC,
               Covered=G$Covered[A,])
   class(out) <- "cccd"
   out
}


dominate.weighted.cccd <- function(G,dxx,alpha=0,igncover=F,verbose=F)
{
   C <- G$Adj
   N <- G$Num
   A <- NULL
   NC <- 0
   nx <- nrow(C)
   weights <- 1/(1+apply(dxx,1,sum))
   if(alpha>=nx) stop("alpha too big: cannot eliminate all the points")
   num <- nx-alpha
   mw <- min(c(weights,0))-1

   covered <- rep(0,nx)

   while(NC<num){
      i <- biggestplace(weights)
	  A <- c(A,i)
	  covered <- covered | C[i,]
	  NC <- sum(covered)
	  #weights[C[i,]>0] <- mw
	  dd <- dxx[,!covered]
	  if(is.matrix(dd))
		 weights <- 1/(1+apply(dd,1,sum))
	  else
		 weights <- 1/(1+dd)
	  weights[A] <- mw
	  if(igncover) weights[covered] <- mw
	  if(verbose) cat("\n",NC)
   }
   nw <- length(A)
   out <- list(Indices=A,Radii=G$R[A],AdjacencyMatrix=G$Adj[A,A],NumCovered=apply(G$Adj[A,],1,sum),
               Covered=G$Covered[A,])
   class(out) <- "cccd"
   out
}

depth.local <- function(dx,R)
{
   nx <- nrow(dx)
   inside.ball <- t(t(dx)<R)
   d1 <- apply(dx*inside.ball,1,sum)
   d2 <- apply(inside.ball,1,sum)
   d3 <- 1/(d1/d2+1)
   d3[d2==1] <- 0
   d3
}

