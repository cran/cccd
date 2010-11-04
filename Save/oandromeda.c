/* Description:
 *
 */

#include <stdio.h>
#include <math.h>

#define MIN(x,y) (((x)<(y))?(x):(y))

/* sorts biggest to smallest */
int cmp(const void *x, const void *y)
{
   double **xx,**yy;

   xx = (double **)x;
   yy = (double **)y;

   if((*xx)[0]>(*yy)[0]) return(-1);
   if((*xx)[0]<(*yy)[0]) return(1);
   return(0);
}

/* sorts smallest to biggest */
int cmp1(const void *x, const void *y)
{
   double *xx,*yy;

   xx = (double *)x;
   yy = (double *)y;

   if((*xx)<(*yy)) return(-1);
   if((*xx)>(*yy)) return(1);
   return(0);
}

/* sorts smallest to biggest */
int cmp2(const void *x, const void *y)
{
   double **xx,**yy;

   xx = (double **)x;
   yy = (double **)y;

   if((*xx)[0]<(*yy)[0]) return(-1);
   if((*xx)[0]>(*yy)[0]) return(1);
   return(0);
}

void andromeda(double *DX,    /* distances between x and x */
               double *DY,    /* distances between y and x */
			   int *NX,       /* number of x points */
			   int *NY,       /* number of y points */
			   int *A,        /* index for witness elements */
			   double *R,     /* radiuses */
			   int *NW,       /* number of witness elements */
			   int *VERBOSE)
{
   int i,j,k,l;
   int nx=*NX;
   int ny=*NY;
   double **dx;
   double **dy;
   double **ord;
   int ncovered;
   int *covered;
   int verbose=*VERBOSE;

   /* 
    * dx is nx by nx
	* dy is ny by nx
	*/
   dx = (double **)calloc(nx,sizeof(double *));
   for(i=0;i<nx;i++){
      dx[i] = DX+i*nx;
   }

   dy = (double **)calloc(ny,sizeof(double *));
   for(i=0;i<ny;i++){
      dy[i] = DY+i*nx; /* note: this IS nx here! */
   }
   /*
   for(j=0;j<ny;j++){
      fprintf(stderr,"%lf ",dy[j][0]);
   }
   */

   covered = (int *)calloc(nx,sizeof(int));

   ord = (double **)calloc(nx,sizeof(double *));
   for(i=0;i<nx;i++){
      ord[i] = (double *)calloc(3,sizeof(double));
	  ord[i][1] = i;
	  ord[i][0] = dy[0][i];
	  for(j=1;j<ny;j++){
	     ord[i][0] = MIN(ord[i][0],dy[j][i]);
	  }
   /* fprintf(stderr,"\n%lf",ord[i][0]); */
   }
   qsort(ord,nx,sizeof(double *),cmp);

   /* fprintf(stderr,"\n%lf %lf %lf %lf\n",ord[0][0],ord[1][0],ord[2][0],ord[3][0]); */

   ncovered=0;

   i = 0; /* number of witness elts so far */
   j = 0; /* position in the ord vector */
   while(ncovered < nx){
      A[i] = (int)ord[j][1];
	  R[i] = ord[j][0];
	  if(verbose) fprintf(stderr,"\nA[%d] = %d, R = %g",i,A[i],R[i]);
	  /* eliminate covered points from consideration */
	  /* first mark all points covered */
	  for(k=0;k<nx;k++){
	     if(dx[k][A[i]]<=R[i]){
			if(covered[k] == 0) ncovered++;
			covered[k] = 1;
		 }
	  }
	  /* then mark points in the ord matrix as covered */
	  for(k=j;k<nx;k++){
	     if(covered[(int)ord[k][1]]){
		    ord[k][2] = 1;
		 }
	  }
	  i++;
	  /* find next uncovered point */
	  if(ncovered<nx){
	     while(ord[j][2]) j++;
	  }
	  if(verbose) {
		 fprintf(stderr,"\n%d of %d covered so far",ncovered,nx);
		 fprintf(stderr,"\nNext vector: %d",j);
	  }
   }
   *NW = i;

   free(dx);
   free(dy);
   for(i=0;i<nx;i++){
	  free(ord[i]);
   }
   free(ord);
   free(covered);
}

void andromedaC(double *DX,    /* distances between x and x */
               double *DY,    /* distances between y and x */
			   int *NX,       /* number of x points */
			   int *NY,       /* number of y points */
			   int *A,        /* index for witness elements */
			   double *R,     /* radiuses */
			   int *C,        /* points covered */
			   int *NW,       /* number of witness elements */
			   int *TIES,     /* break ties by taking largest radius */
			   int *COVER,    /* eliminate covered centers from consideration */
			   int *VERBOSE)
{
   int i,j,k,l;
   int nx=*NX;
   int ny=*NY;
   double **dx;
   double **dy;
   double *radii;
   int ncovered;
   int **covered;
   int *done;
   int *nc;
   int verbose=*VERBOSE;
   int num;
   int maxn,ind;
   int ties=*TIES;
   int cover=*COVER;

   /* 
    * dx is nx by nx
	* dy is ny by nx
	*/
   dx = (double **)calloc(nx,sizeof(double *));
   for(i=0;i<nx;i++){
      dx[i] = DX+i*nx;
   }

   dy = (double **)calloc(ny,sizeof(double *));
   for(i=0;i<ny;i++){
      dy[i] = DY+i*nx; /* note: this IS nx here! */
   }

   radii = (double *)calloc(nx,sizeof(double));
   for(i=0;i<nx;i++){
	  radii[i] = dy[0][i];
	  for(j=1;j<ny;j++){
	     radii[i] = MIN(radii[i],dy[j][i]);
	  }
   }

   /* find the points covered by each ball */
   done = (int *)calloc(nx,sizeof(int));
   nc = (int *)calloc(nx,sizeof(int));
   covered = (int **)calloc(nx,sizeof(int *));
   for(i=0;i<nx;i++){
      covered[i] = (int *)calloc(nx,sizeof(int));
	  for(j=0;j<nx;j++){
		 if(dx[i][j]<=radii[i]){
			*(C+i*nx+j) = covered[i][j] = 1;
			nc[i]++;
		 }
	  }
   }

   ncovered=0;
   i = 0; /* number of witness elts so far */
   /* find the ball covering the (next) most of the uncovered points */
   while(ncovered < nx){

	  /* the if statements for ties and cover could be moved outside of
	   * these loops for a slight efficiency gain at the expense of 
	   * maintaining two copies of each loop, which is fraught with danger.
	   * so I didn't.
	   */
	  maxn = 0;
      for(j=0;j<nx;j++){
	     if(!done[j]){
			/* take the ball that covers the most uncovered points */
			if(nc[j]>maxn){
			   maxn=nc[j];
			   ind=j;
			}
			else if(ties && (nc[j] == maxn)){
			/* handle ties by taking the dude with the largest radius */
			   if(radii[j]>radii[ind]){
				  maxn=nc[j];
				  ind=j;
			   }
			}
		 }
	  }
	  done[ind] = -1;
	  ncovered += maxn;
	  R[i] = radii[ind];
	  A[i] = ind;
	  if(verbose){
	     fprintf(stderr,"\n%d: %d, %d covered (%d), r = %g",
		                    i,A[i],ncovered,maxn,R[i]);
	  }
	  for(k=0;k<nx;k++){
		 if(!done[k]){
			for(l=0;l<nx;l++){
			   if(covered[ind][l]) {
				  /* remove the newley covered points */
				  nc[k] -= covered[k][l];
				  covered[k][l]=0;
				  if(cover) done[l] = -1;
			   }
			}
		 }
	  }
	  i++;
   }
   *NW = i;

   free(dx);
   free(dy);
   free(radii);
   free(done);
   for(i=0;i<nx;i++){
      free(covered[i]);
   }
   free(covered);
}

/* 
 * DX is nx by nx
 * DY is nx by ny
 * note:
 * this is the TRANSPOSE of the matrix in the other versions.
 */
void andromedaCab(double *DX,    /* distances between x and x */
               double *DY,    /* distances between x and y */
			   int *NX,       /* number of x points */
			   int *NY,       /* number of y points */
			   int *A,        /* index for witness elements */
			   double *R,     /* radiuses */
			   int *C,        /* points covered */
			   int *NW,       /* number of witness elements */
			   int *TIES,     /* break ties by taking largest radius */
			   int *COVER,    /* eliminate covered centers from consideration */
			   double *ALPHA, /* tolerance: measure of x not covered */
			   double *BETA,  /* measure of y covered */
			   int *VERBOSE)
{
   int i,j,k,l;
   int nx=*NX;
   int ny=*NY;
   double **dx;
   double *dy;
   int **covered;
   int *done;
   double *radii;
   double *p;
   int k1,k2;
   int ncovered;
   int *nc;
   int verbose=*VERBOSE;
   int num;
   int maxn,ind;
   int ties=*TIES;
   int cover=*COVER;


   k1 = (int)floor((*ALPHA)*nx);
   k2 = (int)floor((*BETA)*ny);

   if(verbose){
      fprintf(stderr,"\nalpha = %g,\tk1 = %d",*ALPHA,k1);
      fprintf(stderr,"\nbeta  = %g,\tk2 = %d",*BETA,k2);
   }

   radii = (double *)calloc(nx,sizeof(double));
   nc = (int *)calloc(nx,sizeof(int));
   done = (int *)calloc(nx,sizeof(int));
   dy = (double *)calloc(ny,sizeof(double));
   p = DY;
   for(i=0;i<nx;i++){
	  memcpy(dy,p,ny*sizeof(double));
	  qsort(dy,ny,sizeof(double),cmp1);
	  radii[i] = dy[k2];
	  p += ny;
   }
   free(dy);
   if(verbose){
      fprintf(stderr,"\nRadii computed:");
	  for(i=0;i<nx;i++){
	     fprintf(stderr,"\n%g",radii[i]);
	  }
   }

   covered = (int **)calloc(nx,sizeof(int *));

   dx = (double **)calloc(nx,sizeof(double *));
   for(i=0;i<nx;i++){
      dx[i] = (double *)calloc(2,sizeof(double));
   }

   /* find the points covered by each ball */
   for(i=0;i<nx;i++){
      covered[i] = (int *)calloc(nx,sizeof(int));
	  for(j=0;j<nx;j++){
		 p = DX+j*nx;
		 for(k=0;k<nx;k++){
			dx[k][0] = *p++;
			dx[k][1] = k;
		 }
		 qsort(dx,nx,sizeof(double *),cmp2);
		 num=0;
	     for(k=0;k<k1+1;k++){
			/* if one of the neighbors of xj is within the ball of xi */
		    if( *(DX+i*nx+(int)dx[k][1]) <= radii[i]) {
			   *(C+i*nx+j) = covered[i][j] = 1;
			   num++;
			}
		 }
		 if(num>0) nc[i]++;
	  }
   }
   for(i=0;i<nx;i++){
      free(dx[i]);
   }
   free(dx);

   if(verbose){
      fprintf(stderr,"\nCover computed:");
	  for(i=0;i<nx;i++){
	     fprintf(stderr,"\n%g %d",radii[i],nc[i]);
	  }
   }


   ncovered=0;
   i = 0; /* number of witness elts so far */
   /* find the ball covering the (next) most of the uncovered points */
   while(ncovered < nx){

	  maxn = 0;
      for(j=0;j<nx;j++){
	     if(!done[j]){
			/* take the ball that covers the most uncovered points */
			if(nc[j]>maxn){
			   maxn=nc[j];
			   ind=j;
			}
			else if(ties && (nc[j] == maxn)){
			/* handle ties by taking the dude with the largest radius */
			   if(radii[j]>radii[ind]){
				  maxn=nc[j];
				  ind=j;
			   }
			}
		 }
	  }
	  done[ind] = -1;
	  ncovered += maxn;
	  R[i] = radii[ind];
	  A[i] = ind;
	  if(verbose){
	     fprintf(stderr,"\n%d: %d, %d covered (%d), r = %g",
		                    i,A[i],ncovered,maxn,R[i]);
	  }
	  for(k=0;k<nx;k++){
		 if(!done[k]){
			for(l=0;l<nx;l++){
			   if(covered[ind][l]) {
				  /* remove the newley covered points */
				  nc[k] -= covered[k][l];
				  covered[k][l]=0;
				  if(cover) done[l] = -1;
			   }
			}
		 }
	  }
	  i++;
   }
   *NW = i;

   free(nc);
   free(done);
   free(radii);
   for(i=0;i<nx;i++){
      free(covered[i]);
   }
   free(covered);
}

/* 
 * D is n by n
 * nc classes
 * P is n by nc
 *
 * Like andromedaCab except for marked observations:
 *
 * Balls cover at most beta of other classes (sum(P[,other_classes])).
 * They can cover a point if covering one of the neighbors
 * (up to sum(P[,obs_class])).
 * Except for the markings on the points, there is no distinction
 * between classes.
 */
void andromedaCabP(double *D, /* distances between observations */
               double *P,     /* Probabilities (markings) for the observations */
			   int *N,        /* number of x points */
			   int *NC,       /* number of classes */
			   double *ALPHA, /* tolerance: measure of x not covered */
			   double *BETA,  /* measure of y covered */
			   int *A,        /* index for witness elements (returned) */
			   double *R,     /* radiuses (returned) */
			   int *CLASS,    /* class assigned to the ball */
			   int *C,        /* points covered (returned) */
			   int *NW,       /* number of witness elements (returned) */
			   int *VERBOSE)
{
   int i,j,k,l;
   int n=*N;                    /* number of observations */
   int nc=*NC;                  /* number of classes */
   int **covered;               /* the A matrix in the paper: points covered */
   double **dists;              /* D as an nxn matrix */
   double **probs;              /* P as an nxnc matrix */
   double alpha=*ALPHA;
   double beta=*BETA;
   double *p,*d;                /* temporary pointers */
   int *c;
   double *pb;                  /* masses for different classes (size is nc) */
   double ***radii;             /* temporary storage for radii, sorted */
   double **balls;              /* balls defined by a radius, a center and a class */
   double *nclass;              /* ``number'' of observations in each class */
   int *num_covered;            /* number covered by ball i */
   int ncovered;                /* total number of points covered (temporary) */
   double r;
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
   num_covered = (int *)calloc(n,sizeof(int));

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
	  balls[i] = (double *)calloc(3,sizeof(double));
      covered[i] = (int *)calloc(n,sizeof(int));
	  for(j=0;j<nc;j++){
	     nclass[j] += probs[i][j];
	  }
   }
   if(verbose){
      fprintf(stderr,"\nMemory allocated");
   }
   /* balls is an array consisting of:
    *
    * balls[][0] = radius
	* balls[][1] = index to center
	* balls[][2] = class associated with the ball
	*/

   /* find the points covered by each ball:
	*
	* 1. For each x, pick the largest ball such that:
	*    after asigning the ball the class with the largest mass
	*    within the ball, the proportion of the mass of all the other
	*    classes combined is no more than beta.
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
		  memset(pb,0,nc*sizeof(double));

		  r = radii[i][ind][0];         /* current radius */

		  l=0;
		  for(k=0;k<n;k++){
		  /* keep track of mass for each class: pbj.
		   * only tally for those observations within the ball.
		   */
			 if(dists[i][k]<r){
				l++;
				for(j=0;j<nc;j++){
				   pb[j] += probs[(int)radii[i][k][1]][j];
				}
			 }
		  }

		  /* find the class for the ball */
		  maxv = pb[0];
		  maxclass = 0;
		  for(j=1;j<nc;j++){
			 if(pb[j]>maxv){
				maxv = pb[j];
				maxclass=j;
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
			 if(i==16) {
				fprintf(stderr,"\n%d %g %d",i,r,maxclass);
				for(k=0;k<nc;k++) fprintf(stderr,"\n%g",probs[i][k]);
			 }
			 /* balls always cover themselves */
			 covered[i][i] = 1;          
			 num_covered[i]++;

			 /* we will break out of the while loop since we found a ball */
			 done=-1;

			 /* now we have a ball for observation xi.
			  * determine which points it covers, using 2. above.
			  */

			 tot = alpha*nclass[maxclass];

			 /* check each point to determine if it is covered by the ball.
			  * only check points which have not already been determined to
			  * be covered.
			  */

			 for(j=0;j<n;j++){
				if(covered[i][j]==0){
				   /* find k = # nearest neighbors to check */
				   k=0;
				   mass=probs[j][maxclass];
				   while(mass<tot*probs[j][maxclass]){
					  k++;
				      mass += probs[(int)radii[j][n-k-1][1]][maxclass];
				   }
				   if(k>=0){
				      for(l=n-1;l>=n-k-1;l--){
					     if(dists[i][(int)radii[j][l][1]] < r){
						    covered[i][j] = 1;
							num_covered[i]++;
							break;
						 }
					  }
				   }
/*
fprintf(stderr,"\n\ti=%d j=%d k=%d mass=%g n=%d,ri=%g,dij=%g",i,j,k,mass,num_covered[i],r,dists[i][j]);
*/
				}
			 }
		  }
		  ind++;
		  if(ind==n) done=-1;
	   }
   }
   c = C;
   for(i=0;i<n;i++){
       memcpy(c,covered[i],n*sizeof(int));
	   c += n;
   }

   if(verbose){
      fprintf(stderr,"\nCover computed\n");
	  /*
	  for(i=0;i<n;i++) fprintf(stderr,"%d ",num_covered[i]);
	  */
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
		 if(num_covered[j]>maxn){
			maxn=num_covered[j];
			ind=j;
		 }
		 else if(num_covered[j] == maxn){
		 /* handle ties by taking the ball with the largest radius */
			if(balls[j][0]>balls[ind][0]){
			   maxn=num_covered[j];
			   ind=j;
			}
		 }
	  }
	  if(maxn<=0){
	     fprintf(stderr,"\nImpossible: no points to cover and we haven't stopped\n");
		 exit(7);
	  }
	  ncovered += maxn;
	  R[i] = balls[ind][0];
	  A[i] = ind;
	  CLASS[i] = (int)balls[ind][2];
	  balls[ind][0] = 0.0;
	  if(verbose){
	     fprintf(stderr,"\n%d: %d, %d covered (%d), r = %g, class = %d",
		                    i,A[i]+1,ncovered,maxn,R[i],CLASS[i]);
	  }
	  for(k=0;k<n;k++){
		 for(l=0;l<n;l++){
			if(covered[ind][l]) {
			   /* remove the newley covered points */
			   num_covered[k] -= covered[k][l];
			   covered[k][l]=0;
			}
		 }
	  }
	  i++;
	  if(verbose){
	     fprintf(stderr,"\n%d points covered so far",ncovered);
	  }
   }
   *NW = i;
   for(i=0;i<n;i++){
	  free(balls[i]);
	  free(covered[i]);
   }
   free(balls);
   free(covered);
   free(num_covered);
}
