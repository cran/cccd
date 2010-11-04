#include <stdio.h>
#include <math.h>

#define MIN(x,y) (((x)<(y))?(x):(y))

/* sorts biggest to smallest */
extern int cmp(const void *x, const void *y);

/* sorts smallest to biggest */
extern int cmp1(const void *x, const void *y);

/* sorts smallest to biggest */
extern int cmp2(const void *x, const void *y);

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
			   int *SHRINK,   /* shrink the balls to reduce over coverage of Y */
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
   int shrink=*SHRINK;
   double maxxd;
   double temp;


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
	  if(shrink){
		 maxxd = 0.0;
		 /* find the farthest x covered by the ball */
		 for(j=0;j<nx;j++){
			temp = *(DX+i*nx+j);
			if(temp<radii[i]){
			   maxxd = (maxxd>temp)?maxxd:temp;
			}
		 }
		 /* shrink the ball down so that it covers only those x's */
		 for(j=0;j<k2;j++){
		    if(dy[j]>maxxd) {
			   maxxd = dy[j];
			   break;
			}
		 }
		 radii[i] = maxxd;
	  }
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

