
/* This software developed by David Marchette
 * Naval Surface Warfare Center, Dahlgren
 * Division, Code B10.  It may be freely used and 
 * modified. All warranties, express or implied, are
 * disclaimed. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <values.h>

#define MIN(x,y) (((x)<(y))?(x):(y))

/* sorts biggest to smallest */
extern int cmp(const void *x, const void *y);

/* sorts smallest to biggest */
extern int cmp1(const void *x, const void *y);

/* sorts smallest to biggest */
extern int cmp2(const void *x, const void *y);

/* 
 * DX is nx by nx
 * DY is nx by ny
 *
 * Returns 
 * radii in R,
 * covered points in C.
 */
void cccd(double *DX,            /* distances between x and x */
               double *DY,       /* distances between x and y */
			   int *NX,          /* number of x points */
			   int *NY,          /* number of y points */
			   double *R,        /* radiuses (returned) */
			   int *C,           /* points covered (returned) */
			   int *K,           /* measure of y covered, AKA Beta */
			   int *SHRINK,      /* shrink the radii back to first covered x */
			   int *VERBOSE)
{
   int i,j,k,l;
   int nx=*NX;
   int ny=*NY;
   double **dx;
   double *dy;
   int *covered;
   double *p;
   int verbose=*VERBOSE;
   int num;
   int maxn,ind;
   double maxxd;
   double temp;
   int kk;
   int shrink=*SHRINK;


   k = *K;

   dy = (double *)calloc(ny,sizeof(double));
   p = DY;
   for(i=0;i<nx;i++){
	  covered = C+i*nx;
	  memset(covered,0,sizeof(int)*nx);
	  memcpy(dy,p,ny*sizeof(double));
	  qsort(dy,ny,sizeof(double),cmp1);
	  if(k>=ny) {
		 R[i] = MAXFLOAT;
		 fprintf(stderr,"\n%g",R[i]);
	  }
	  else{
		 R[i] = dy[k];
	  }
	  maxxd = dy[0];
	  /* find the farthest x covered by the ball */
	  for(j=0;j<nx;j++){
		 temp = *(DX+i*nx+j);
		 if(temp<R[i]){
			maxxd = (maxxd>temp)?maxxd:temp;
			covered[j] = 1;
		 }
	  }
	  /* shrink the ball down so that it covers only those x's */
	  if(shrink){
		 for(j=0;j<k+1;j++){
			if(dy[j]>=maxxd) {
			   maxxd = dy[j];
			   break;
			}
		 }
		 /* shrink the ball so it covers only those x's */
		 R[i] = maxxd;
	  }
	  p += ny;
   }
   free(dy);
}

