dominate.greedy <- function(g,weight=NULL)
{
   od <- degree(g,mode="out")+1
   S <- NULL
   A <- get.adjacency(g)
   diag(A) <- 0
   n <- nrow(A)
   covered <- rep(0,n)
   while(sum(covered)<n){
      i <- which.max(od)
		if(!is.null(weight)){
		   qq <- weight
			qq[od != od[i]] <- -Inf
			i <- which.max(qq)
		}
      covered[A[i,]>0] <- 1
		covered[i] <- 1
      S <- c(S,i)
      A[,covered>0] <- 0
		h <- graph.adjacency(A,mode="directed")
      od <- degree(h,mode="out")+1-covered
   }
   S
}

dominate.byR <- function(g)
{
   dominate.greedy(g,weight=g$R)
}

dominate.random.sample <- function(g)
{
	n <- vcount(g)
   dominate.greedy(g,weight=runif(n))
}

dominate <- function(g,method="greedy")
{
   METHODS = c("greedy","sample","byradius")
   method <- pmatch(tolower(method),METHODS)
   if(is.na(method)){
      stop("invalid method")
   }
	S <- NA
   if(method==1) S <- dominate.greedy(g)
   else if(method==2) S <- dominate.random.sample(g)
	else if(method==3) S <- dominate.byR(g)
   S
}

