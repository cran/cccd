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
   S-1
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

dominate <- function(g,method="greedy",verbose=TRUE)
{
   METHODS = c("greedy","sample","R","systematic")
   method <- pmatch(method,METHODS)
   if(is.na(method)){
      stop("invalid method")
   }
	S <- NA
   if(method==1) S <- dominate.greedy(g)
   else if(method==2) S <- dominate.random.sample(g)
	else if(method==3) S <- dominate.byR(g)
	else if(method==4) S <- dominate.systematic(g,verbose=verbose)
   S
}

dominate.systematic <- function(g,S,verbose=TRUE)
{
	x <- degree(g,mode="out")
   A <- get.adjacency(g)
   n <- nrow(A)
	if(max(x) == nrow(A)) return(which.max(x))
   if(missing(S)) S <- dominate.greedy(g)
   gamma <- length(S)
   if(verbose){
      cat("Starting with gamma =",gamma,"\n")
   }
   if(gamma <= 2) return(S)
   zeros <- degree(g,mode="in")==0
   Z <- which(zeros)-1
   W <- which(zeros==FALSE)-1
   if((gamma-1-length(Z))==0){
      if(length(unique(unlist(neighborhood(g,order=1,nodes=Z))))==n){
         return(Z)
      }
      else{
         return(S)
      }
   }
   if(sum(zeros)<n){
      numchoices <- choose(length(W),gamma-1-length(Z))
      if(verbose){
         cat("This may require","(",length(W),",",gamma-1-length(Z),") =",
             numchoices,"steps\n")
      }
      m <- gamma-1-length(Z)
      txt <- paste("for(i1 in 1:",length(W)-(m-1),")",sep="")
      arg <- paste("W[i1]",sep="")
      j <- 1
      fnd <- 0
		S1 <- NULL
      if(m>1){
         for(i in (m-2):0){
            txt <- paste(txt," if(!fnd) for(i",j+1," in (i",j,"+1):",length(W)-i,")",sep="")
            arg <- paste(arg,",W[","i",j+1,"]",sep="")
            j <- j+1
         }
      }
      txt <- paste(txt,
         "{",
         "if(!fnd){",
         " S1 <- c(Z,",arg,");",
			"nbd <- unique(unlist(neighborhood(g,order=1,nodes=S1)));",
         "if(length(nbd)==",n,"){fnd <- 1; break }",
         "}",
         "}",
      sep="")
      eval(parse(text=txt))
      if(fnd){
         S <- S1
         if(verbose){
            cat("Recursing:",S,"\n")
         }
         Q <- Recall(g,S)
			if(length(Q)<length(S)) S <- Q
      }
   }
   S
}
