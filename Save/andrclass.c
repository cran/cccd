#include <stdio.h>
#include <math.h>

/* sorts biggest to smallest */
extern int cmp(const void *x, const void *y);

/* 
 * D is n by n distance matrix
 * nc classes
 * P is n by nc probabilities (markings) for the observation classes
 *
 * Like andromedaCab except for marked observations:
 *
 * Balls cover at most beta of other classes (sum(P[,other_classes])).
 * They can cover a point if covering one of the neighbors
 * (up to sum(P[,obs_class])).
 * Except for the markings on the points, there is no distinction
 * between classes.
 */
void andromedaCabP(double *D, /* distances between observations                   */
               double *P,     /* Probabilities (markings) for the observations    */
			   int *N,        /* number of x points                               */
			   int *NC,       /* number of classes                                */
			   double *ALPHA, /* tolerance: measure of x not covered              */
			   double *BETA,  /* measure of y covered                             */
			   int *BYCENTER, /* if true, use center of ball for class assignment */
			   int *A,        /* index for witness elements (returned) 			  */
			   double *R,     /* radiuses (returned) 			 			   	  */
			   int *CLASS,    /* class assigned to the ball 			 		  */
			   int *C,        /* points covered (returned if not null) 			  */
			   int *NW,       /* number of witness elements (returned) 			  */
			   int *VERBOSE)
{
   int i,j,k,l;
   int n=*N;                    /* number of observations                    */
   int nc=*NC;                  /* number of classes                         */
   int **covered;               /* the A matrix in the paper: points covered */
   double **dists;              /* D as an nxn matrix                        */
   double **probs;              /* P as an nxnc matrix                       */
   double alpha=*ALPHA;
   double beta=*BETA;
   double *p,*d;                /* temporary pointers                        */
   int *c;
   double *pb;                  /* masses for different classes (size is nc) */
   double ***radii;             /* temporary storage for radii, sorted       */
   double **balls;              /* balls defined by a radius, a center and a class */
   double *nclass;              /* ``number'' of observations in each class   */
   int ncovered;                /* total number of points covered (temporary) */
   double r;
   int bycenter=*BYCENTER;
   int verbose=*VERBOSE;
   int num;
   int ind;
   int done;
   int maxclass;
   double maxv;
   double mass;
   double tot;
   int maxn;


   if(verbose){
      fprintf(stderr,"\nalpha = %g",*ALPHA);
      fprintf(stderr,"\nbeta  = %g",*BETA);
   }

   pb = (double *)calloc(nc,sizeof(double));
   radii = (double ***)calloc(n,sizeof(double **));
   balls = (double **)calloc(n,sizeof(double *));

   nclass = (double *)calloc(nc,sizeof(double));

   covered = (int **)calloc(n,sizeof(int *));

   dists = (double **)calloc(n,sizeof(double *));
   probs = (double **)calloc(n,sizeof(double *));
   d = D;
   p = P;
   for(i=0;i<n;i++){
      dists[i] = d;
	  probs[i] = p;
	  d += n;
	  p += nc;
	  radii[i] = (double **)calloc(n,sizeof(double *));
	  for(j=0;j<n;j++){
	     radii[i][j] = (double *)calloc(2,sizeof(double));
		 radii[i][j][0] = dists[i][j];
		 radii[i][j][1] = j;
	  }
	  qsort(radii[i],n,sizeof(double *),cmp);
	  balls[i] = (double *)calloc(4,sizeof(double));
      covered[i] = (int *)calloc(n,sizeof(int));
	  for(j=0;j<nc;j++){
	     nclass[j] += probs[i][j];
	  }
   }
   if(verbose){
      fprintf(stderr,"\nMemory allocated\n");
   }
   /* balls is an array consisting of:
    *
    * balls[][0] = radius
	* balls[][1] = index to center
	* balls[][2] = class associated with the ball
	* balls[][3] = number of points covered
	*/

   /* find the points covered by each ball:
	*
	* 1. For each x, pick the largest ball such that:
	*    after asigning the ball the class with the largest P
	*    the proportion of the mass of all the other
	*    classes combined for all the points in the ball is no more than beta.
	* 2. Given a ball j, and an observation i, we determine if the ball
	*    ``covers'' i by the following rule:
	*    Let k be the largest integer such that:
	*    sum(probs[i][l]) < alpha*Nc*P[i][c]
	*    where the sum is from 0 to k, and Nc is the sum of probs[][c]
	*    for the class assigned to the ball.
	*    Then i is covered by the ball iff one of its k nearest neighbors
	*    is covered:  dist[j][l] < rj.
	*    If P[i][c] = 1, then we are saying that we'll allow neighbors up
	*    to 100alpha percent of the total. For P[i][c] < 1, we are using
	*    the ``amount'' of class c ``in'' the observation to adjust this.
	*    So if the observation has P[i][c] = 0 (is definitely not a class
	*    c observation), then we only consider the observation itself (k=0).
    */
   for(i=0;i<n;i++){
	  /* start at the biggest ball, and check condition 1. above.
	   * if it fails, move to a smaller ball. continue until either
	   * pass, or reach the end. information goes in the variable balls.
	   */
	   done = 0;
	   ind = 0;                   /* index into radii */

	   while(!done){
		  /* clear the pb array */
		  memset(pb,0,nc*sizeof(double));

		  r = radii[i][ind][0];         /* current radius */

		  for(k=ind+1;k<n;k++){
		  /* keep track of mass for each class: pbj.
		   * only tally for those observations within the ball.
		   */
			 for(j=0;j<nc;j++){
				pb[j] += probs[(int)radii[i][k][1]][j];
			 }
		  }

		  /* find the class for the ball */
		  if(bycenter){
			 maxv = probs[i][0];
			 maxclass = 0;
			 for(j=1;j<nc;j++){
				if(probs[i][j]>maxv){
				   maxv = probs[i][j];
				   maxclass=j;
				}
			 }
		  }
		  else {
			 maxv = pb[0];
			 maxclass = 0;
			 for(j=1;j<nc;j++){
				if(pb[j]>maxv){
				   maxv = pb[j];
				   maxclass=j;
				}
			 }
		  }

		  /* find the total mass for all the other classes */
		  mass =0.0;
		  tot = 0.0;
		  for(j=0;j<nc;j++){
		     mass += pb[j];
			 tot += nclass[j];
		  }
		  mass -= pb[maxclass];
		  tot -= nclass[maxclass];

		  /* now check the condition: 
		   * mass of all the other classes is less than
		   * beta*Nc, that is, 100*beta percent of the mass.
		   */
		  if(mass <= beta*tot){
		     balls[i][0] = r;
			 balls[i][1] = i;
			 balls[i][2] = maxclass;

			 /*
			 fprintf(stderr,"\nball %d, class %d, radius %g",i,maxclass,r);
			 */

			 /* we will break out of the while loop since we found a ball */
			 done=-1;

			 /*
			 if(maxclass == 0){
			    fprintf(stderr,"\nr=%g, i = %d",r,i);
			 }
			 */

			 /* now we have a ball for observation xi.
			  * determine which points it covers, using 2. above.
			  */

			 tot = alpha*nclass[maxclass];

			 /* check each point to determine if it is covered by the ball. */

			 for(j=0;j<n;j++){
				/* find k = # nearest neighbors to check */
				k=0;
				mass=probs[j][maxclass];
				
				while(mass<tot*probs[j][maxclass]){
				   k++;
				   mass += probs[(int)radii[j][n-k-1][1]][maxclass];
				}
				for(l=n-1;l>=n-k-1;l--){
				   if(dists[i][(int)radii[j][l][1]] < r){
					  covered[i][j] = 1;
					  balls[i][3]++;
					  break;
				   }
				}
/*
fprintf(stderr,"\n\ti=%d j=%d k=%d mass=%g n=%d,ri=%g,dij=%g",i,j,k,mass,balls[i][3],r,dists[i][j]);
*/
			 }
			 /*
			 fprintf(stderr,", covered %d",(int)balls[i][3]);
			 */
		  }
		  ind++;
		  if(ind==n) done=-1;
	   }
	   if(verbose){
	      if((i+1)%10==0) fprintf(stderr," %d",i+1);
	   }
   }
   if(C){
	  c = C;
	  for(i=0;i<n;i++){
		  memcpy(c,covered[i],n*sizeof(int));
		  c += n;
	  }
   }

   if(verbose){
      double maxr,minr;
	  int numzero,numnear;
      fprintf(stderr,"\nCover computed\nRadii:\n");
	  minr = maxr = balls[0][0];
	  numzero=0;
	  numnear=0;
	  for(i=1;i<n;i++){
		 minr = (minr>balls[i][0])?balls[i][0]:minr;
		 maxr = (maxr<balls[i][0])?balls[i][0]:maxr;
		 if(balls[i][0]==0.0) numzero++;
		 if(balls[i][0]<1E-8) numnear++;
	  }
	  fprintf(stderr,"Min Radius:  %g  Max Radius: %g\n",minr,maxr);
	  fprintf(stderr,"Number 0 radii:  %d  Number < 1E-8: %d\n",numzero,numnear);
   }

   free(pb);
   free(dists);
   free(probs);
   free(nclass);
   for(i=0;i<n;i++){
      for(j=0;j<n;j++){
	     free(radii[i][j]);
	  }
	  free(radii[i]);
   }
   free(radii);

   if(verbose){
      fprintf(stderr,"\nFreed memory (radii)");
   }


   ncovered=0;
   i = 0; /* number of witness elts so far */
   /* find the ball covering the (next) most of the uncovered points */
   while(ncovered < n){
      
	  maxn = 0;
	  ind = 0;
      for(j=0;j<n;j++){
		 /* take the ball that covers the most uncovered points */
		 if(balls[j][3]>maxn){
			maxn=(int)balls[j][3];
			ind=j;
		 }
		 else if((int)balls[j][3] == maxn){
		 /* handle ties by taking the ball with the largest radius */
			if(balls[j][0]>balls[ind][0]){
			   ind=j;
			}
		 }
	  }
	  if(maxn<=0){
	     fprintf(stderr,"\nShouldn't happen: no points to cover and we haven't stopped.");
	     fprintf(stderr,"\nCould be caused by coincident observations of more than one class.\n");
		 goto endgame;
	  }
	  ncovered += maxn;
	  R[i] = balls[ind][0];
	  A[i] = ind;
	  CLASS[i] = (int)balls[ind][2];

	  /* eliminate the ball from further consideration */
	  balls[ind][0] = -1.0;             
	  /* balls[ind][3] = 0.0;             */

	  if(verbose){
	     fprintf(stderr,"\n%d: %d, %d covered (%d), r = %g, class = %d",
		                    i,A[i]+1,ncovered,maxn,R[i],CLASS[i]);
	  }
	  for(l=0;l<n;l++){
		 if(covered[ind][l]) {
			for(k=0;k<n;k++){
			   /* remove the newley covered points */
			   if(covered[k][l]){
				  balls[k][3]--;
				  covered[k][l]=0;
			   }
			}
		 }
	  }
	  i++;
	  if(verbose){
	     fprintf(stderr,"\n%d points covered so far",ncovered);
	  }
   }
   endgame:
   if(verbose){
      fprintf(stderr,"\n");
   }
   *NW = i;
   for(i=0;i<n;i++){
	  free(balls[i]);
	  free(covered[i]);
   }
   free(balls);
   free(covered);
}
