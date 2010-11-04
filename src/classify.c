
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

#define PRECISION 1E-8

void classify_andromeda_class(
               double *D,     /* distances between observations and witness elts: nXnw       */
			   int *N,        /* number of x observations                                    */
			   int *CLASS,    /* class assigned to the witness elts: vector of length nw     */
			   int *NC,       /* number of classes */
			   double *R,     /* radii: vector of length nw                                  */
			   int *NW,       /* number of witness elements                                  */
			   int *out,      /* output class assignments (-1 means no call)                 */
			   int *DOUBLES,  /* output flag for observations in more than one class' balls  */
			   double *MAXR)  /* if d/r>maxr no call is made                                 */
{
   int i,j;
   int n=*N;                    /* number of observations */
   int nw=*NW;                  /* number of witnesses    */
   double **dists;              /* D as an nxn matrix     */
   double *d;                   /* temporary pointer      */
   int *di;                     /* temporary pointer      */
   double maxr=*MAXR;
   int nc=*NC;
   int ind;
   int **doubles;
   double minv;


   /* doubles is now nXnc */
   if(DOUBLES) {
	  doubles = (int **)calloc(n,sizeof(double *));
	  di = DOUBLES;
	  for(i=0;i<n;i++){
	     doubles[i] = di;
		 di += nc;
		 memset(doubles[i],0,n*sizeof(int));
	  }
   }
   d = D;
   for(i=0;i<n;i++){
	  minv = d[0]/R[0];
	  ind = 0;
	  for(j=1;j<nw;j++){
	     if(minv>d[j]/R[j]){
		    minv = d[j]/R[j];
			ind = j;
		 }
	  }
	  if(minv<maxr){
		 out[i] = CLASS[ind];
	  }
	  else{
	     out[i] = -1;
	  }
	  if(doubles){
	     for(j=0;j<nw;j++){
			if(CLASS[j] != CLASS[ind]){
			   if(d[j]<R[j]){
				  if((R[j]-d[j])>PRECISION){ /* kludge for round-off errors */
					 doubles[i][CLASS[j]]++;
				  }
			   }
			}
		 }
	  }
	  d += nw;
   }
   if(DOUBLES) free(doubles);
}
#undef PRECISION
